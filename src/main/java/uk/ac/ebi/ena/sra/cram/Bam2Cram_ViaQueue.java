package uk.ac.ebi.ena.sra.cram;

import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.Random;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.zip.GZIPOutputStream;

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

import uk.ac.ebi.ena.sra.cram.bam.BAMFileQueryQueues;
import uk.ac.ebi.ena.sra.cram.bam.SAMRecordIteratorJob;
import uk.ac.ebi.ena.sra.cram.bam.Sam2CramRecordFactory;
import uk.ac.ebi.ena.sra.cram.format.CramCompression;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;
import uk.ac.ebi.ena.sra.cram.impl.ByteArraySequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.impl.CramHeaderIO;
import uk.ac.ebi.ena.sra.cram.impl.SequentialCramWriter;
import uk.ac.ebi.ena.sra.cram.spot.PairedTemplateAssembler;
import uk.ac.ebi.ena.sra.cram.stats.CramStats;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.FileConverter;

public class Bam2Cram_ViaQueue {

	private static Logger log;

	public static void main(String[] args) throws Exception {
		log = Logger.getLogger(Bam2Cram_ViaQueue.class);

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
		long dumpRecords = 0;
		int maxReadLength = Integer.MAX_VALUE;
		double coverageModifier = 1.0;
		boolean skipPerfectMatch = false;

		log.info("Input BAM file: " + bamFile);
		SAMFileReader
				.setDefaultValidationStringency(ValidationStringency.SILENT);
		SAMFileReader samFileReader = new SAMFileReader(bamFile);

		FileOutputStream cramFOS = null;
		OutputStream bos = null;
		ByteArrayOutputStream cramBAOS = null;
		boolean wrapInGzip = params.gzip;
		if (outputCramFile != null) {
			log.info("Output file: " + outputCramFile.getAbsolutePath());
			cramFOS = new FileOutputStream(outputCramFile);
			if (wrapInGzip)
				bos = new BufferedOutputStream(new GZIPOutputStream(cramFOS));
			else
				bos = new BufferedOutputStream(cramFOS);
		} else {
			cramBAOS = new ByteArrayOutputStream();
			if (wrapInGzip)
				bos = new BufferedOutputStream(new GZIPOutputStream(cramBAOS));
			else
				bos = new BufferedOutputStream(cramBAOS);
		}
		if (params.roundTripCheck)
			log.info("Round trip check is on, expect slight delays.");

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
		cramHeader.setVersion("0.3");
		cramHeader
				.setReferenceSequences(new ArrayList<CramReferenceSequence>());

		for (SAMSequenceRecord samRF : header.getSequenceDictionary()
				.getSequences()) {
			cramHeader.getReferenceSequences().add(
					new CramReferenceSequence(samRF.getSequenceName(), samRF
							.getSequenceLength()));
		}

		CramHeaderIO.write(cramHeader, bos);

		ExecutorService es = Executors.newFixedThreadPool(1);

		for (String seqName : seqNames) {

			ReferenceSequence sequence = referenceSequenceFile
					.getSequence(seqName);
			byte[] refBases = referenceSequenceFile.getSubsequenceAt(
					sequence.getName(), 1, sequence.length()).getBases();
			Utils.capitaliseAndCheckBases(refBases);

			byte[] refStart = new byte[50];
			System.arraycopy(refBases, 0, refStart, 0, refStart.length);
			log.info("Reference sequence " + seqName + ", starts with "
					+ new String(refStart));

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
			CramStats stats = new CramStats();

			counter = 1L;
			long nofSoftClips = 0L;
			long softClipLength = 0L;
			long samRecordCounter = 0L;
			Sam2CramRecordFactory cramRecordFactory = new Sam2CramRecordFactory(
					refBases);
			PairedTemplateAssembler assembler = new PairedTemplateAssembler(
					params.spotAssemblyAlignmentHorizon,
					params.spotAssemblyRecordsHorizon);

			final SAMRecord stopSAMRecord = SAMRecordIteratorJob.STOP_SAMRECORD;
			// final SAMRecordIterator iterator =
			// samFileReader.queryOverlapping(
			// seqName, 0, 0);
			// if (!iterator.hasNext())
			// continue;

			BAMFileQueryQueues qq = new BAMFileQueryQueues(es, bamFile);
			Queue<SAMRecord> samRecordQueue = qq.getQueueForQuery(seqName, 0,
					0, true);
			// new Thread(new Runnable() {
			//
			// @Override
			// public void run() {
			// while (iterator.hasNext()) {
			// try {
			// samRecordQueue.put(iterator.next());
			// } catch (InterruptedException e) {
			// // TODO Auto-generated catch block
			// e.printStackTrace();
			// break;
			// }
			// }
			// samRecordQueue.add(stopSAMRecord);
			// }
			// }).start();

			long nullCounter = 0;
			SAMRecord samRecord = null;
			try {
				while ((samRecord = samRecordQueue.poll()) != stopSAMRecord) {
					if (samRecord == null) {
						// Thread.sleep(1000) ;
						nullCounter++;
						continue;
					}
					if (skipRecords > 0 && samRecordCounter++ < skipRecords)
						continue;
					// if (samRecord.getReadUnmappedFlag())
					// continue;

					assembler.addSAMRecord(samRecord);

					if ((samRecord = assembler.nextSAMRecord()) == null)
						continue;

					int distance = assembler.distanceToNextFragment();

					if (samRecord.getReadLength() > maxReadLength)
						Utils.changeReadLength(samRecord, maxReadLength);

					for (CigarElement ce : samRecord.getCigar()
							.getCigarElements()) {
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
					cramRecord.setReadMapped(true);
					if (distance > 0) {
						cramRecord.setLastFragment(false);
						cramRecord.setRecordsToNextFragment(distance);
					} else
						cramRecord.setLastFragment(true);
					// if (skipPerfectMatch && cramRecord.isPerfectMatch())
					// continue;
					// cramRecord.setLastFragment(true);
					stats.addRecord(cramRecord);
					if (counter++ >= maxRecords)
						break;
				}
			} finally {
				System.out.println("null counter=" + nullCounter);
				qq.stop(samRecordQueue);
			}
			log.info("Nof soft clips: " + nofSoftClips);
			log.info("Nof of soft clipped bases: " + softClipLength);

			stats.adjustBlock(block);
			// iterator.close();

			log.info(block);
			writer.write(block);
			writer.flush();
			if (cramBAOS != null)
				log.info("Block size: " + cramBAOS.size());

			counter = 1;
			samRecordCounter = 0L;
			Random random = new Random();

			// final SAMRecordIterator iterator2 =
			// samFileReader.queryOverlapping(
			// seqName, 0, 0);
			samRecordQueue.clear();
			samRecordQueue = qq.getQueueForQuery(seqName, 0, 0, true);
			// new Thread(new Runnable() {
			//
			// @Override
			// public void run() {
			// samRecordQueue.clear();
			// while (iterator2.hasNext()) {
			// try {
			// samRecordQueue.put(iterator2.next());
			// } catch (InterruptedException e) {
			// // TODO Auto-generated catch block
			// e.printStackTrace();
			// break;
			// }
			// }
			// samRecordQueue.add(stopSAMRecord);
			// }
			// }).start();

			assembler.clear();
			nullCounter = 0;
			try {
				while ((samRecord = samRecordQueue.poll()) != stopSAMRecord) {
					// while (iterator2.hasNext()) {
					// samRecord = iterator2.next() ;

					if (samRecord == null) {
						// Thread.sleep(1000) ;
						nullCounter++;
						continue;
					}

					CramRecord cramRecord;
					if (skipRecords > 0 && samRecordCounter++ < skipRecords)
						continue;

					if (coverageModifier < 1.0
							&& random.nextFloat() > coverageModifier)
						continue;

					assembler.addSAMRecord(samRecord);

					if ((samRecord = assembler.nextSAMRecord()) == null)
						continue;
					int distance = assembler.distanceToNextFragment();

					// if (samRecord.getReadUnmappedFlag())
					// continue;

					if (samRecord.getReadLength() > maxReadLength)
						Utils.changeReadLength(samRecord, maxReadLength);

					cramRecord = cramRecordFactory.createCramRecord(samRecord);
					cramRecord.setReadMapped(true);
					if (distance > 0) {
						cramRecord.setLastFragment(false);
						cramRecord.setRecordsToNextFragment(distance);
					} else
						cramRecord.setLastFragment(true);
					// if (skipPerfectMatch && cramRecord.isPerfectMatch())
					// continue;
					// cramRecord.setLastFragment(true);
					try {
						if (counter < dumpRecords
								|| (counter > failingRecord - 5 && counter < failingRecord + 5)) {
							byte[] ref = new byte[50];
							System.arraycopy(refBases,
									samRecord.getAlignmentStart() - 1, ref, 0,
									ref.length);
							log.info(new String(ref));
							log.info(counter + ": " + cramRecord);
						}
					} catch (Throwable e1) {
						e1.printStackTrace();
						break;
					}
					bases += cramRecord.getReadLength();

					try {
						if (params.roundTripCheck)
							writer.writeAndCheck(cramRecord);
						else
							writer.write(cramRecord);
					} catch (Exception e) {
						System.err.println(cramRecord);
						e.printStackTrace();
						throw e;
					}
					if (counter++ >= maxRecords)
						break;
				}
			} finally {
				qq.stop(samRecordQueue);
				System.out.println("null counter=" + nullCounter);
			}
			totalCounter += counter;
			// iterator.close();
			writer.flush();
			if (log.isDebugEnabled())
				writer.dump();
		}

		es.shutdownNow();

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

		@Parameter(names = { "--reference-fasta-file" }, converter = FileConverter.class, required = true)
		File referenceFasta;

		@Parameter(names = { "--output-cram-file" }, converter = FileConverter.class)
		File outputFile = null;

		@Parameter(names = { "--max-records" })
		long maxRecords = Long.MAX_VALUE;

		@Parameter
		List<String> sequences;

		@Parameter(names = { "-h", "--help" }, description = "Print help and quit")
		boolean help = false;

		@Parameter(names = { "--round-trip-check" })
		boolean roundTripCheck = false;

		@Parameter(names = { "--gzip" })
		boolean gzip = false;

		@Parameter(names = { "--record-horizon" })
		int spotAssemblyRecordsHorizon = 10000;

		@Parameter(names = { "--alignment-horizon" })
		int spotAssemblyAlignmentHorizon = 10000;
	}
}
