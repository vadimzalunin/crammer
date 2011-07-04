package uk.ac.ebi.ena.sra.cram;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.LinkedBlockingQueue;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;
import uk.ac.ebi.ena.sra.cram.impl.ByteArraySequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.impl.CramHeaderIO;
import uk.ac.ebi.ena.sra.cram.impl.ReadFeatures2Cigar;
import uk.ac.ebi.ena.sra.cram.impl.RestoreBases;
import uk.ac.ebi.ena.sra.cram.impl.RestoreQualityScores;
import uk.ac.ebi.ena.sra.cram.impl.SequentialCramReader;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.FileConverter;

public class Cram2Bam {
	private static Logger log;

	private static void setupLogger() {
		// PatternLayout layout = new PatternLayout(
		// "%d{ABSOLUTE} %5p %c{1}:%L - %m%n");
		//
		// ConsoleAppender appender = new ConsoleAppender(layout, "System.err");
		// appender.setThreshold(Level.INFO);
		//
		// Logger.getRootLogger().addAppender(appender);
		// Logger.getRootLogger().setLevel(Level.ALL);
		log = Logger.getLogger(Cram2Bam.class);
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

		ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory
				.getReferenceSequenceFile(params.reference);
		FileInputStream cramFIS = new FileInputStream(params.cramFile);

		cram2bam(referenceSequenceFile, cramFIS, params.outputFile,
				params.maxRecords);

	}

	public static void cram2bam(ReferenceSequenceFile referenceSequenceFile,
			InputStream cramInputStream, File outputBamFile, long maxRecords)
			throws Exception {

		DataInputStream cramDIS = new DataInputStream(cramInputStream);
		CramHeader cramHeader = CramHeaderIO.read(cramDIS);

		BAMFileWriter writer = new BAMFileWriter(outputBamFile);
		writer.setSortOrder(SortOrder.coordinate, true);
		SAMFileHeader header = new SAMFileHeader();
//		SAMFileWriter out = new SAMFileWriterFactory().makeSAMWriter(header,
//				true, System.out);
		for (CramReferenceSequence cs : cramHeader.getReferenceSequences()) {
			SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord(
					cs.getName(), cs.getLength());
			header.addSequence(samSequenceRecord);
		}
		writer.setHeader(header);

		long time1 = System.currentTimeMillis();
		NEXT_BLOCK: while (true) {
			SequentialCramReader reader = new SequentialCramReader(cramDIS,
					null);
			CramRecordBlock readBlock = reader.readBlock();
			if (readBlock == null)
				break NEXT_BLOCK;
			log.info(readBlock);
			SAMSequenceRecord sequence = header.getSequence(readBlock
					.getSequenceName());
			if (sequence == null) {
				sequence = new SAMSequenceRecord(readBlock.getSequenceName(),
						readBlock.getSequenceLength());
				header.addSequence(sequence);
			}

			String seqName = readBlock.getSequenceName();

			ReferenceSequence nextSequence = referenceSequenceFile
					.getSequence(seqName);
			byte[] refBases = referenceSequenceFile.getSubsequenceAt(
					nextSequence.getName(), 1, nextSequence.length())
					.getBases();
			Utils.capitaliseAndCheckBases(refBases);
			ByteArraySequenceBaseProvider provider = new ByteArraySequenceBaseProvider(
					refBases);
			reader.setReferenceBaseProvider(provider);

			int counter = 1;
			long start = 0L;
			CramRecord cramRecord = null;
			RestoreBases restoreBases = new RestoreBases();
			restoreBases.setProvider(provider);
			restoreBases.setSequenceName(seqName);
			RestoreQualityScores restoreQualityScores = new RestoreQualityScores();
			ReadFeatures2Cigar readFeatures2Cigar = new ReadFeatures2Cigar();

			ArrayList<CramRecord> records = new ArrayList<CramRecord>(
					(int) readBlock.getRecordCount());
			readBlock.setRecords(records);
			for (long i = 0; i < readBlock.getRecordCount(); i++) {
				try {
					cramRecord = reader.readRecord();
					records.add(cramRecord);

					if (counter++ >= maxRecords)
						break;

				} catch (Exception e) {
					log.error("Failed to read record: " + i);
					log.error(cramRecord.toString());
					log.error(e);
					throw e;
				}
			}

			counter = 1;
			Map<Integer, Integer> indexes = new TreeMap<Integer, Integer>();
			for (CramRecord record : records) {
				cramRecord = record;
				start += cramRecord.getAlignmentStart();
				if (!cramRecord.isLastFragment()) {
					if (cramRecord.getRecordsToNextFragment() + counter > records
							.size())
						cramRecord.setLastFragment(true);
					else {
						indexes.put((int) (counter + cramRecord
								.getRecordsToNextFragment()), counter);

					}
				}

				SAMRecord samRecord = new SAMRecord(header);

				Integer index = indexes.remove(counter);
				if (index != null) {
					samRecord.setReadName(String.valueOf(index.intValue())
							+ ".2");
				} else {
					if (cramRecord.isLastFragment())
						samRecord.setReadName(String.valueOf(counter));
					else {
						samRecord.setReadName(String.valueOf(counter) + ".1");
					}
				}

				samRecord.setAlignmentStart((int) cramRecord
						.getAlignmentStart());
				if (cramRecord.isReadMapped()) {
					samRecord.setReadBases(restoreBases
							.restoreReadBases(cramRecord));
					samRecord.setBaseQualities(restoreQualityScores
							.restoreQualityScores(cramRecord));
				} else {
					samRecord.setReadBases(cramRecord.getReadBases());
					samRecord.setBaseQualities(cramRecord.getQualityScores());
				}
				samRecord.setCigar(readFeatures2Cigar.getCigar(
						cramRecord.getReadFeatures(),
						(int) cramRecord.getReadLength()));
				samRecord.setReadUnmappedFlag(false);
				samRecord.setReferenceName(readBlock.getSequenceName());
				writer.addAlignment(samRecord);
//				out.addAlignment(samRecord);
				if (counter++ >= maxRecords)
					break NEXT_BLOCK;
			}
			log.info(start);
		}

//		out.close();
		writer.close();
		long time2 = System.currentTimeMillis();
		log.info("Decoded in: " + (time2 - time1) + " millis");
	}

	@Parameters(commandDescription = "CRAM to BAM conversion")
	static class Params {
		@Parameter(names = { "--input-cram-file" }, converter = FileConverter.class)
		File cramFile;

		@Parameter(names = { "--max-sequences" })
		int maxSequences = Integer.MAX_VALUE;

		@Parameter(names = { "--max-records" })
		long maxRecords = Long.MAX_VALUE;

		@Parameter(names = { "--reference-fasta-file" }, required = true, converter = FileConverter.class)
		File reference;

		@Parameter(names = { "--output-bam-file" }, converter = FileConverter.class)
		File outputFile;

		@Parameter(names = { "-h", "--help" }, description = "Print help and quit")
		boolean help = false;
	}
}
