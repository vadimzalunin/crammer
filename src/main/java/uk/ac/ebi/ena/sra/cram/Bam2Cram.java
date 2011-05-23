package uk.ac.ebi.ena.sra.cram;

import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.bam.Sam2CramRecordFactory;
import uk.ac.ebi.ena.sra.cram.format.CramCompression;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;
import uk.ac.ebi.ena.sra.cram.impl.ByteArraySequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.impl.CramHeaderIO;
import uk.ac.ebi.ena.sra.cram.impl.SequentialCramWriter;
import uk.ac.ebi.ena.sra.cram.stats.CramStats;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.FileConverter;

public class Bam2Cram {

	private static Logger log;

	public static void main(String[] args) throws Exception {
		log = Logger.getLogger(Bam2Cram.class) ;
		
		Params params = new Params();
		JCommander jc = new JCommander(params);
		jc.parse(args);

		File bamFile = params.bamFile;
		File refFile = params.referenceFasta;
		File outputCramFile = params.outputFile;

		long maxRecords = Long.MAX_VALUE ;
		List<String> seqNames = params.sequences;
		long dumpRecords = 0;
		int maxReadLength = Integer.MAX_VALUE;
		double coverageModifier = 1.0;
		boolean skipPerfectMatch = false;

		log.info("Input BAM file: " + bamFile);
		SAMFileReader
				.setDefaultValidationStringency(ValidationStringency.SILENT);
		SAMFileReader samFileReader = new SAMFileReader(bamFile);

		FileOutputStream cramFOS = null;
		BufferedOutputStream bos = null;
		ByteArrayOutputStream cramBAOS = null;
		if (outputCramFile != null) {
			log.info("Output file: " + outputCramFile.getAbsolutePath());
			cramFOS = new FileOutputStream(outputCramFile);
			bos = new BufferedOutputStream(cramFOS);
		} else {
			cramBAOS = new ByteArrayOutputStream();
			bos = new BufferedOutputStream(cramBAOS);
		}

		SAMFileHeader header = samFileReader.getFileHeader();

		if (seqNames == null)
			seqNames = new ArrayList<String>();
		if (seqNames.isEmpty()) {
			int maxSequences = Integer.MAX_VALUE;
			for (SAMSequenceRecord seq : header.getSequenceDictionary()
					.getSequences()) {
				seqNames.add(seq.getSequenceName());
				if (seqNames.size() >= maxSequences)
					break;
			}
		}

		ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory
				.getReferenceSequenceFile(refFile);

		long counter = 1L;
		long totalCounter = 1L;
		long bases = 0L;
		long time1 = System.currentTimeMillis();

		long failingRecord = Long.MAX_VALUE;
		long skipRecords = 0;
		CramHeader cramHeader = new CramHeader();
		cramHeader.setVersion("0.2");
		cramHeader
				.setReferenceSequences(new ArrayList<CramReferenceSequence>());

		for (SAMSequenceRecord samRF : header.getSequenceDictionary()
				.getSequences()) {
			cramHeader.getReferenceSequences().add(
					new CramReferenceSequence(samRF.getSequenceName(), samRF
							.getSequenceLength()));
		}

		CramHeaderIO.write(cramHeader, bos);

		for (String seqName : seqNames) {

			ReferenceSequence sequence = referenceSequenceFile
					.getSequence(seqName);
			byte[] refBases = referenceSequenceFile.getSubsequenceAt(
					sequence.getName(), 1, sequence.length()).getBases();
			byte[] refStart = new byte[50];
			System.arraycopy(refBases, 0, refStart, 0, refStart.length);
			log.info("Reference sequence " + seqName + ", starts with "
					+ new String(refStart));

			for (int i = 0; i < refBases.length; i++) {
				switch (refBases[i]) {
				case 'A':
				case 'C':
				case 'G':
				case 'T':
				case 'N':
					break;

				default:
					throw new RuntimeException("Illegal base at " + i + ": "
							+ refBases[i]);
				}
			}

			ByteArraySequenceBaseProvider provider = new ByteArraySequenceBaseProvider(
					refBases);
			SequentialCramWriter writer = new SequentialCramWriter(bos,
					provider);

			CramRecordBlock block = new CramRecordBlock();
			block.setPositiveStrandBasePositionReversed(false);
			block.setNegativeStrandBasePositionReversed(true);
			block.setSequenceName(seqName);
			block.setSequenceLength(refBases.length);
			block.setCompression(new CramCompression());

			samFileReader = new SAMFileReader(bamFile);
			SAMRecordIterator recordIterator = samFileReader.queryOverlapping(
					seqName, 0, 0);
			CramStats stats = new CramStats();

			if (!recordIterator.hasNext())
				continue;

			counter = 1L;
			long nofSoftClips = 0L;
			long softClipLength = 0L;
			long samRecordCounter = 0L;
			Sam2CramRecordFactory cramRecordFactory = new Sam2CramRecordFactory(
					refBases);
			while (recordIterator.hasNext()) {
				SAMRecord samRecord = recordIterator.next();
				if (skipRecords > 0 && samRecordCounter++ < skipRecords)
					continue;
				if (samRecord.getReadUnmappedFlag())
					continue;

				if (samRecord.getReadLength() > maxReadLength)
					Utils.changeReadLength(samRecord, maxReadLength);

				for (CigarElement ce : samRecord.getCigar().getCigarElements()) {
					switch (ce.getOperator()) {
					case S:
						nofSoftClips++;
						softClipLength += ce.getLength();
						break;

					default:
						break;
					}
				}

				CramRecord cramRecord = cramRecordFactory
						.createCramRecord(samRecord);
				if (skipPerfectMatch && cramRecord.isPerfectMatch())
					continue;
				cramRecord.setLastFragment(true);
				stats.addRecord(cramRecord);
				if (counter++ >= maxRecords)
					break;
			}
			log.info("Nof soft clips: " + nofSoftClips);
			log.info("Nof of soft clipped bases: " + softClipLength);

			stats.adjustBlock(block);
			recordIterator.close();

			log.info(block);
			writer.write(block);
			writer.flush();
			if (cramBAOS != null)
				log.info("Block size: " + cramBAOS.size());

			counter = 1;
			samRecordCounter = 0L;
			SAMRecordIterator iterator = samFileReader.queryOverlapping(
					seqName, 0, 0);
			Random random = new Random();
			while (iterator.hasNext()) {
				SAMRecord samRecord = iterator.next();

				CramRecord cramRecord;
				if (skipRecords > 0 && samRecordCounter++ < skipRecords)
					continue;
				if (samRecord.getReadUnmappedFlag())
					continue;

				if (coverageModifier < 1.0
						&& random.nextFloat() > coverageModifier)
					continue;

				if (samRecord.getReadLength() > maxReadLength)
					Utils.changeReadLength(samRecord, maxReadLength);

				cramRecord = cramRecordFactory.createCramRecord(samRecord);
				if (skipPerfectMatch && cramRecord.isPerfectMatch())
					continue;
				cramRecord.setLastFragment(true);
				try {
					if (counter < dumpRecords
							|| (counter > failingRecord - 5 && counter < failingRecord + 5)) {
						byte[] ref = new byte[50];
						System.arraycopy(refBases,
								samRecord.getAlignmentStart() - 1, ref, 0,
								ref.length);
						log.info(new String(ref));
						log.info(counter + ": " + cramRecord);
						log.info(new String(Utils.restoreBases(cramRecord,
								provider, seqName)));
					}
				} catch (Throwable e1) {
					e1.printStackTrace();
					break;
				}
				bases += cramRecord.getReadLength();

				try {
					writer.write(cramRecord);
				} catch (Exception e) {
					System.err.println(cramRecord);
					e.printStackTrace();
					throw e;
				}
				if (counter++ >= maxRecords)
					break;
			}
			totalCounter += counter;
			iterator.close();
			writer.flush();
		}

		long time2 = System.currentTimeMillis();
		bos.close();

		log.info("Written " + totalCounter + " reads");
		long cramLength = cramBAOS == null ? outputCramFile.length() : cramBAOS
				.size();
		log.info("File size: " + cramLength);
		log.info("In: " + (time2 - time1) + " millis");
		log.info("Bytes per read: " + (float) cramLength / totalCounter);
		log.info("Bits per base: " + (float) cramLength * 8 / bases);
	}

	@Parameters(commandDescription = "BAM to CRAM converter.")
	static class Params {
		@Parameter(names = { "--input-bam-file" }, converter = FileConverter.class, required = true)
		File bamFile;

		@Parameter(names = { "--reference-fasta" }, converter = FileConverter.class, required = true)
		File referenceFasta;

		@Parameter(names = { "--output-cram-file" }, converter = FileConverter.class)
		File outputFile = null;

		@Parameter
		List<String> sequences;
	}
}
