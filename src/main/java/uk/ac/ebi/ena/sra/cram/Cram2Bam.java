package uk.ac.ebi.ena.sra.cram;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;

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
	private static Logger log = Logger.getLogger(Cram2Bam.class);

	private static InputStream createCramInputStream(File file)
			throws IOException {
		FileInputStream fis = new FileInputStream(file);

		// gzip magic:
		if (fis.read() == 31 && fis.read() == 139)
			return new GZIPInputStream(new BufferedInputStream(
					new FileInputStream(file)));

		return new BufferedInputStream(new FileInputStream(file));
	}

	public static void main(String[] args) throws Exception {
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
		InputStream cramIS = createCramInputStream(params.cramFile);

		convert(referenceSequenceFile, cramIS, params.outputFile,
				params.maxRecords);

	}

	public static void convert(ReferenceSequenceFile referenceSequenceFile,
			InputStream cramInputStream, File outputBamFile, long maxRecords)
			throws Exception {

		DataInputStream cramDIS = new DataInputStream(cramInputStream);
		CramHeader cramHeader = CramHeaderIO.read(cramDIS);

		BAMFileWriter writer = new BAMFileWriter(outputBamFile);
		writer.setSortOrder(SortOrder.coordinate, true);
		SAMFileHeader header = new SAMFileHeader();

		OutputStream samOS = new FlushNotCloseOutputStream(System.out);
		SAMFileWriter out = new SAMFileWriterFactory().makeSAMWriter(header,
				true, samOS);
		for (CramReferenceSequence cs : cramHeader.getReferenceSequences()) {
			SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord(
					cs.getName(), cs.getLength());
			header.addSequence(samSequenceRecord);
		}
		writer.setHeader(header);

		long indexOfFirstRecordInTheBlock = 0L;
		CramRecordBlock prevBlock = null;
		long prevAlStart = 0L;

		long time1 = System.currentTimeMillis();

		NEXT_BLOCK: while (true) {
			SequentialCramReader reader = new SequentialCramReader(cramDIS,
					null);
			CramRecordBlock readBlock = reader.readBlock();
			if (readBlock == null)
				break NEXT_BLOCK;

			if (prevBlock != null) {
				if (!prevBlock.getSequenceName().equals(
						readBlock.getSequenceName())) {
					prevAlStart = 0L;
				}
			}

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
			Utils.capitaliseAndCheckBases(refBases, false);
			ByteArraySequenceBaseProvider provider = new ByteArraySequenceBaseProvider(
					refBases);
			reader.setReferenceBaseProvider(provider);

			int counter = 1;
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
					if (cramRecord != null)
						log.error(cramRecord.toString());
					log.error(e);
					throw e;
				}
			}

			counter = 1;

			Map<Integer, Integer> indexes = new TreeMap<Integer, Integer>();

			for (CramRecord record : records) {
				cramRecord = record;
				if (!cramRecord.isLastFragment()) {
					if (cramRecord.getRecordsToNextFragment() + counter > records
							.size())
						cramRecord.setLastFragment(true);
					else {
						indexes.put((int) (counter + cramRecord
								.getRecordsToNextFragment()), counter);

					}
				}
				Integer index = indexes.remove(counter);

				SAMRecord samRecord = new SAMRecord(header);
				if (index != null) {
					samRecord.setReadName(String
							.valueOf(indexOfFirstRecordInTheBlock
									+ index.intValue())
							+ ".2");
					samRecord.setSecondOfPairFlag(true);
					samRecord.setReadPairedFlag(true) ;
				} else {
					if (cramRecord.isLastFragment()) {
						samRecord
								.setReadName(String
										.valueOf(indexOfFirstRecordInTheBlock
												+ counter));
						samRecord.setReadPairedFlag(false) ;
					}
					else {
						samRecord
								.setReadName(String
										.valueOf(indexOfFirstRecordInTheBlock
												+ counter)
										+ ".1");
						samRecord.setFirstOfPairFlag(true);
						samRecord.setReadPairedFlag(true) ;
					}
				}

				if (cramRecord.isReadMapped()) {
					samRecord.setAlignmentStart((int) cramRecord
							.getAlignmentStart());
					samRecord.setReadBases(restoreBases
							.restoreReadBases(cramRecord));
					byte[] scores = restoreQualityScores
							.restoreQualityScores(cramRecord);
					injectQualityScores(scores, samRecord);
					prevAlStart = samRecord.getAlignmentStart();
				} else {
					samRecord.setAlignmentStart((int) prevAlStart);
					samRecord.setReadBases(cramRecord.getReadBases());
					byte[] scores = cramRecord.getQualityScores();
					injectQualityScores(scores, samRecord);
				}
				samRecord.setCigar(readFeatures2Cigar.getCigar(
						cramRecord.getReadFeatures(),
						(int) cramRecord.getReadLength()));
				samRecord.setReadUnmappedFlag(!cramRecord.isReadMapped());
				samRecord.setReadNegativeStrandFlag(cramRecord
						.isNegativeStrand());
				samRecord.setReferenceName(readBlock.getSequenceName());
				try {
					writer.addAlignment(samRecord);
					// out.addAlignment(samRecord);
				} catch (IllegalArgumentException e) {
					log.error("Offensive block: " + readBlock.toString());
					log.error("Offensive CRAM record: " + cramRecord.toString());
					log.error("Offensive SAM record: " + samRecord.toString());
					log.error("SAM record al start="
							+ samRecord.getAlignmentStart());
					throw e;
				}

				if (counter++ >= maxRecords) {
					indexOfFirstRecordInTheBlock += records.size();
					break NEXT_BLOCK;
				}
			}
			indexOfFirstRecordInTheBlock += records.size();
			prevBlock = readBlock;
		}

		out.close();
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

	private static class FlushNotCloseOutputStream extends OutputStream {
		private OutputStream delegate;

		public FlushNotCloseOutputStream(OutputStream delegate) {
			super();
			this.delegate = delegate;
		}

		@Override
		public void write(int b) throws IOException {
			delegate.write(b);
		}

		@Override
		public void write(byte[] b) throws IOException {
			delegate.write(b);
		}

		@Override
		public void write(byte[] b, int off, int len) throws IOException {
			delegate.write(b, off, len);
		}

		@Override
		public void close() throws IOException {
			delegate.flush();
			delegate = null;
		}

	}

	private static final void injectQualityScores(byte[] scores,
			SAMRecord record) {
		if (scores == null || scores.length == 0) {
			injectNullQualityScores(record);
			return;
		}
		final byte nullQS = -1;
		final byte asciiOffset = 33;
		final byte defaultQS = '?';
		final byte space = 32;

		boolean nonDefaultQsFound = false;
		for (int i = 0; i < scores.length; i++)
			if (scores[i] != space) {
				nonDefaultQsFound = true;
				break;
			}

		if (!nonDefaultQsFound) {
			injectNullQualityScores(record);
			return;
		}

		for (int i = 0; i < scores.length; i++) {
			scores[i] -= asciiOffset;
			if (scores[i] == nullQS)
				scores[i] = defaultQS - asciiOffset;
		}

		record.setBaseQualities(scores);
	}

	private static final void injectNullQualityScores(SAMRecord record) {
		record.setBaseQualities(new byte[0]);
	}
}
