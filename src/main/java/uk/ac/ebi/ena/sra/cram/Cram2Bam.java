package uk.ac.ebi.ena.sra.cram;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;

import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;

import uk.ac.ebi.ena.sra.cram.format.CramFormatException;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;
import uk.ac.ebi.ena.sra.cram.format.compression.CramCompressionException;
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
		PatternLayout layout = new PatternLayout(
				"%d{ABSOLUTE} %5p %c{1}:%L - %m%n");

		ConsoleAppender appender = new ConsoleAppender(layout, "System.err");
		appender.setThreshold(Level.INFO);

		Logger.getRootLogger().addAppender(appender);
		Logger.getRootLogger().setLevel(Level.ALL);
		log = Logger.getLogger(Cram2Bam.class);
	}

	public static void main(String[] args) throws IOException,
			CramFormatException, CramCompressionException {
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
				.getReferenceSequenceFile(new File(
						"y:/Data/psyringae/psyringae.fa"));
		FileInputStream cramFIS = new FileInputStream(
				"y:/Data/psyringae/psyringae.cram");

		cram2bam(referenceSequenceFile, cramFIS, new File(
				"y:/Data/psyringae/cram2bam.bam"));

	}

	public static void cram2bam(ReferenceSequenceFile referenceSequenceFile,
			InputStream cramInputStream, File outputBamFile)
			throws CramCompressionException, IOException {

		DataInputStream cramDIS = new DataInputStream(cramInputStream);
		CramHeader cramHeader = CramHeaderIO.read(cramDIS);
		long maxRecords = Long.MAX_VALUE;

		BAMFileWriter writer = new BAMFileWriter(outputBamFile);
		writer.setSortOrder(SortOrder.coordinate, true);
		SAMFileHeader header = new SAMFileHeader();
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
			Utils.capitaliseAndCheckBases(refBases) ;
			ByteArraySequenceBaseProvider provider = new ByteArraySequenceBaseProvider(
					refBases);
			reader.setReferenceBaseProvider(provider);

			long counter = 1;
			long start = 0L;
			CramRecord cramRecord = null;
			RestoreBases restoreBases = new RestoreBases();
			restoreBases.setProvider(provider);
			restoreBases.setSequenceName(seqName);
			RestoreQualityScores restoreQualityScores = new RestoreQualityScores();
			ReadFeatures2Cigar readFeatures2Cigar = new ReadFeatures2Cigar();
			for (long i = 0; i < readBlock.getRecordCount(); i++) {
				try {
					cramRecord = reader.readRecord();

					start += cramRecord.getAlignmentStart();

					SAMRecord samRecord = new SAMRecord(header);
					samRecord.setAlignmentStart((int) cramRecord
							.getAlignmentStart());
					samRecord.setReadBases(restoreBases
							.restoreReadBases(cramRecord));
					samRecord.setBaseQualities(restoreQualityScores
							.restoreQualityScores(cramRecord));
					samRecord.setReadName(String.valueOf(counter));
					samRecord.setCigar(readFeatures2Cigar.getCigar(
							cramRecord.getReadFeatures(),
							(int) cramRecord.getReadLength()));
					samRecord.setReadUnmappedFlag(false);
					samRecord.setReferenceName(readBlock.getSequenceName()) ;
					writer.addAlignment(samRecord);

					if (counter++ >= maxRecords)
						break NEXT_BLOCK;
				} catch (Exception e) {
					log.error("Failed to read record: " + i);
					log.error(cramRecord.toString());
					log.error(e);
					break NEXT_BLOCK;
				}
			}
			log.info(start);
		}

		writer.close();
		long time2 = System.currentTimeMillis();
		log.info("Decoded in: " + (time2 - time1) + " millis");
	}

	@Parameters(commandDescription = "CRAM to BAM conversion")
	static class Params {
		@Parameter(names = { "--input-cram" }, converter = FileConverter.class)
		File cramFile;

		@Parameter(names = { "--max-sequences" })
		int maxSequences = 0;

		@Parameter(names = { "--max-records-per-sequence" })
		long maxRecordsPerSequence = 0;

		@Parameter(names = { "--reference" }, required = true, converter = FileConverter.class)
		File reference;

		@Parameter(names = { "--output-file" }, converter = FileConverter.class)
		File outputFile;
		
		@Parameter(names = { "-h", "--help" }, description = "Print help and quit")
		boolean help = false;
	}
}
