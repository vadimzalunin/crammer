package uk.ac.ebi.ena.sra.cram;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
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

import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;

import uk.ac.ebi.ena.sra.cram.bam.Sam2CramRecordFactory;
import uk.ac.ebi.ena.sra.cram.format.CramCompression;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;
import uk.ac.ebi.ena.sra.cram.impl.ByteArraySequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.impl.CramHeaderIO;
import uk.ac.ebi.ena.sra.cram.impl.RestoreBases;
import uk.ac.ebi.ena.sra.cram.impl.SequentialCramReader;
import uk.ac.ebi.ena.sra.cram.impl.SequentialCramWriter;
import uk.ac.ebi.ena.sra.cram.stats.CramStats;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.FileConverter;

public class CramWriterReadTest {
	private static Logger log;

	private static void setupLogger() {
		PatternLayout layout = new PatternLayout(
				"%d{ABSOLUTE} %5p %c{1}:%L - %m%n");

		ConsoleAppender appender = new ConsoleAppender(layout, "System.err");
		appender.setThreshold(Level.INFO);

		Logger.getRootLogger().addAppender(appender);
		Logger.getRootLogger().setLevel(Level.ALL);
		log = Logger.getLogger(CramWriterReadTest.class);
	}

	public static void main(String[] args) throws Exception {
		setupLogger();

		Params params = new Params();
		JCommander jc = new JCommander(params);
		jc.parse(args);
		
		if (args.length == 0 || params.help) {
			StringBuilder sb = new StringBuilder();
			sb.append("\n");
			jc.usage(sb);

			System.out.println(sb.toString());
			return;
		}

		File bamFile = params.bamFile;
		File refFile = params.referenceFasta;
		File outputCramFile = params.outputFile;

		long maxRecords = params.maxRecords;
		List<String> seqNames = params.sequences;
		long dumpRecords = params.dumpRecords;
		int maxReadLength = params.maxRecordLength;
		double coverageModifier = params.coverageModifier;
		boolean skipPerfectMatch = false;
		boolean testAllRecords = params.roundTripCheck;

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
		List<CramRecord> readCramRecords = new ArrayList<CramRecord>();
		List<SAMRecord> samRecords = new ArrayList<SAMRecord>();
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
			Utils.capitaliseAndCheckBases(refBases) ;
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
				// if (!samRecord.getReadNegativeStrandFlag())
				// continue;
				// if (!samRecord.getCigarString().contains("I")
				// || !samRecord.getCigarString().contains("D"))
				// continue;

				// if (!samRecord.getCigarString().contains("S")) continue ;

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
				// CramRecord cramRecord = CramRecordStaticFactory.newRecord3(
				// samRecord, refBases);
				if (skipPerfectMatch && cramRecord.isPerfectMatch())
					continue;
				cramRecord.setLastFragment(true);
				if (testAllRecords) {
					readCramRecords.add(cramRecord);
					samRecords.add(samRecord);
				}
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

				// if (!samRecord.getReadNegativeStrandFlag())
				// continue;

				// if (!samRecord.getCigarString().contains("S")) continue ;
				// if (!samRecord.getCigarString().contains("I")
				// || !samRecord.getCigarString().contains("D"))
				// continue;

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
						// log.info(Utils.toString(
						// String.valueOf(counter),
						// Utils.restoreBases(cramRecord, provider,
						// block.getSequenceName()), cramRecord));
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

		counter = 0;

		InputStream cramIS = cramBAOS == null ? new BufferedInputStream(
				new FileInputStream(outputCramFile))
				: new ByteArrayInputStream(cramBAOS.toByteArray());

		CramHeader cramHeader2 = CramHeaderIO.read(cramIS);
		log.info(cramHeader2.toString());

		DataInputStream cramDIS = new DataInputStream(cramIS);

		time1 = System.currentTimeMillis();
		for (String seqName : seqNames) {

			ReferenceSequence nextSequence = referenceSequenceFile
					.getSequence(seqName);
			byte[] refBases = referenceSequenceFile.getSubsequenceAt(
					nextSequence.getName(), 1, nextSequence.length())
					.getBases();
			ByteArraySequenceBaseProvider provider = new ByteArraySequenceBaseProvider(
					refBases);
			SequentialCramReader reader = new SequentialCramReader(cramDIS,
					provider);
			CramRecordBlock readBlock = reader.readBlock();
			log.info(readBlock);

			counter = 1;
			long start = 0L;
			CramRecord cramRecord = null;
			RestoreBases restoreBases = new RestoreBases();
			restoreBases.setProvider(provider);
			restoreBases.setSequenceName(seqName);
			for (long i = 0; i < readBlock.getRecordCount(); i++) {
				try {
					cramRecord = reader.readRecord();
					if (i < dumpRecords
							|| (i > failingRecord - 5 && i < failingRecord + 5)) {
						log.info(i + ": " + cramRecord);
					}

					if (testAllRecords) {
						if (!cramRecord.equals(readCramRecords.get((int) i))) {
							log.error("Expecting record: ");
							log.error(readCramRecords.get((int) i));
							log.error(samRecords.get((int) i).format());
							log.error("Found record: ");
							log.error(cramRecord.toString());
							break;
						}
						byte[] samBases = samRecords.get((int) i)
								.getReadBases();
						byte[] cramBases = restoreBases
								.restoreReadBases(cramRecord);
						if (!Arrays.equals(samBases, cramBases)) {
							log.error("Read bases don't match: ");
							log.error(readCramRecords.get((int) i));
							log.error(new String(restoreBases
									.restoreReadBases(readCramRecords
											.get((int) i))));
							log.error(samRecords.get((int) i).format());
							log.error(new String(samRecords.get((int) i)
									.getReadBases()));
							log.error("Found record: ");
							log.error(cramRecord.toString());
							log.error(new String(restoreBases
									.restoreReadBases(cramRecord)));
						}
					}
					start += cramRecord.getAlignmentStart();
					if (counter++ >= maxRecords)
						break;
				} catch (Exception e) {
					e.printStackTrace();
					System.err.println("Failed to read record: " + i);
					System.err.println(cramRecord.toString());
					break;
				}

			}
			log.info(start);
		}
		time2 = System.currentTimeMillis();
		log.info("Decoded in: " + (time2 - time1) + " millis");
	}

	@Parameters(commandDescription = "BAM to CRAM converter, test and development.")
	static class Params {
		@Parameter(names = { "--input-bam-file" }, converter = FileConverter.class, required = true)
		File bamFile;

		@Parameter(names = { "--max-records" })
		long maxRecords = Long.MAX_VALUE;

		@Parameter(names = { "--max-record-length" })
		int maxRecordLength = Integer.MAX_VALUE;
		
		@Parameter(names = { "--dump-records" })
		int dumpRecords = Integer.MAX_VALUE;
		
		@Parameter(names = { "--coverage-modifier" })
		float coverageModifier = 1.0F;

		@Parameter(names = { "--reference-fasta" }, converter = FileConverter.class, required = true)
		File referenceFasta;

		@Parameter(names = { "--output-cram-file" }, converter = FileConverter.class)
		File outputFile = null;

		@Parameter(names = { "--round-trip-check" })
		boolean roundTripCheck = false;

		@Parameter
		List<String> sequences;
		
		@Parameter(names = { "-h", "--help" }, description = "Print help and quit")
		boolean help = false;
	}
}
