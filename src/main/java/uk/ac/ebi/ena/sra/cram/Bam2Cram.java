package uk.ac.ebi.ena.sra.cram;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
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

import uk.ac.ebi.ena.sra.cram.encoding.JavaBinaryCramCodec;
import uk.ac.ebi.ena.sra.cram.format.CramCompression;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.CramRecordFactory;
import uk.ac.ebi.ena.sra.cram.impl.EncodingConstants;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.FileConverter;

public class Bam2Cram {

	public static void main(String[] args) throws IOException {
		PatternLayout layout = new PatternLayout(
				"%d{ABSOLUTE} %5p %c{1}:%L - %m%n");

		ConsoleAppender appender = new ConsoleAppender(layout, "System.err");
		appender.setThreshold(Level.ALL);

		Logger.getRootLogger().addAppender(appender);
		Logger.getRootLogger().setLevel(Level.ALL);

		Logger log = Logger.getLogger(Bam2Cram.class.getSimpleName());

		log.debug("Starting...");

		Params params = new Params();
		JCommander jc = new JCommander(params);
		jc.parse(args);

		ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory
				.getReferenceSequenceFile(new File(params.reference));

		SAMFileReader
				.setDefaultValidationStringency(ValidationStringency.SILENT);

		SAMFileReader reader = new SAMFileReader(params.bamFile);

		if (!reader.hasIndex()) {
			System.err.println("Index: not found");
			return;
		}

		SAMFileHeader header = reader.getFileHeader();

		List<String> seqNames = new ArrayList<String>();

		DataOutputStream os = null;
		if (params.outputFile != null)
			os = new DataOutputStream(new BufferedOutputStream(
					new FileOutputStream(params.outputFile), 1024 * 1024));

		for (SAMSequenceRecord seq : header.getSequenceDictionary()
				.getSequences()) {
			seqNames.add(seq.getSequenceName());
			if (params.maxSequences > 0
					&& seqNames.size() >= params.maxSequences)
				break;
		}

		long totalRecordCounter = 0;
		CramHeader cramHeader = new CramHeader();
		cramHeader.setBlocks(new ArrayList<CramRecordBlock>());

		log.debug("Free nmemory: " + Runtime.getRuntime().freeMemory());
		for (String seqName : seqNames) {
			ReferenceSequence sequence = referenceSequenceFile
					.getSequence(seqName);
			byte[] refSequence = referenceSequenceFile.getSubsequenceAt(
					seqName, 0, sequence.length()).getBases();
			log.info(seqName + ": " + refSequence.length + " bases.");

			CramRecordBlock block = new CramRecordBlock();
			CramCompression compression = new CramCompression();
			compression.setInReadPosGolombLogM(params.intraReadGolombLogM);
			compression.setInSeqPosGolombLogM(params.interReadGolombLogM);
			compression.setDelLengthGolombLogM(params.delLengthGolombLogM);
			block.setCompression(compression);

			List<CramRecord> blockRecords = new ArrayList<CramRecord>();
			block.setRecords(blockRecords);
			int seqRecordCounter = 0;
			SAMRecordIterator iterator;
			iterator = reader.queryOverlapping(seqName, 0, 0);

			SAMRecord record;

			long factoryTime = 0;
			while (iterator.hasNext()) {
				if (seqRecordCounter > params.maxRecordsPerSequence
						&& params.maxRecordsPerSequence > 0)
					break;

				record = iterator.next();
				if (record.getReadUnmappedFlag())
					continue;

				long time = System.currentTimeMillis();
				CramRecord cramRecord = CramRecordFactory.newRecord3(record,
						refSequence);
				factoryTime += System.currentTimeMillis() - time;

				if (params.skipPerfectMatch && cramRecord.isPerfectMatch())
					continue;
				blockRecords.add(cramRecord);
				block.setReadLength(record.getReadLength());
				seqRecordCounter++;
			}
			totalRecordCounter += seqRecordCounter;
			cramHeader.getBlocks().add(block);
			block.setRecordCount(blockRecords.size());
			log.debug("CRAM records creation time: " + factoryTime);
		}

		log.debug("Free nmemory: " + Runtime.getRuntime().freeMemory());
		log.debug("Saving...");
		JavaBinaryCramCodec codec = new JavaBinaryCramCodec(null);
		codec.write(os, cramHeader);
		log.debug("Done.");
	}

	@Parameters(commandDescription = "BAM to CRAM converter.")
	static class Params {
		@Parameter(names = { "--bam-file" }, converter = FileConverter.class, required = true)
		File bamFile;

		@Parameter(names = { "--max-sequences" })
		int maxSequences = 0;

		@Parameter(names = { "--max-records-per-sequence" })
		long maxRecordsPerSequence = 0;

		@Parameter(names = { "--reference" }, required = true)
		String reference;

		@Parameter(names = { "--output-file" }, converter = FileConverter.class)
		File outputFile;

		@Parameter(names = { "--skip-perfect-match" })
		boolean skipPerfectMatch = false;

		@Parameter(names = { "--inter-read-golomb-logm" })
		int interReadGolombLogM = EncodingConstants.INTER_READ_GOLOMB_LOG2M;

		@Parameter(names = { "--intra-read-golomb-logm" })
		int intraReadGolombLogM = EncodingConstants.INTRA_READ_GOLOMB_LOG2M;

		@Parameter(names = { "--del-length-golomb-logm" })
		int delLengthGolombLogM = EncodingConstants.DELETION_LENGTH_GOLOMB_LOG2M;
	}
}
