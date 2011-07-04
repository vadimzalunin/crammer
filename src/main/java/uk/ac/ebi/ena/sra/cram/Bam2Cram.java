package uk.ac.ebi.ena.sra.cram;

import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.bam.Sam2CramRecordFactory;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;
import uk.ac.ebi.ena.sra.cram.impl.CramWriter;
import uk.ac.ebi.ena.sra.cram.spot.PairedTemplateAssembler;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.FileConverter;

public class Bam2Cram {

	private static Logger log = Logger.getLogger(Bam2Cram.class);

	private File bamFile;
	private File refFile;
	private SAMFileReader samReader;
	private CramWriter cramWriter;
	private PairedTemplateAssembler assembler;
	private Sam2CramRecordFactory cramRecordFactory;
	private List<CramReferenceSequence> sequences;
	private OutputStream os;
	private SequenceBaseProvider provider;

	private final boolean roundTripCheck;
	private final int spotAssemblyAlignmentHorizon;
	private final int spotAssemblyRecordsHorizon;

	private ReferenceSequenceFile referenceSequenceFile;

	private final long maxRecords;
	private long recordCount;
	private long unmappedRecordCount;
	private long baseCount;

	private final int maxBlockSize;
	private long recordsBitLength;

	public Bam2Cram(File bamFile, OutputStream os, File refFile,
			long maxRecords, boolean roundTripCheck,
			int spotAssemblyAlignmentHorizon, int spotAssemblyRecordsHorizon,
			int maxBlockSize) {
		this.bamFile = bamFile;
		this.os = os;
		this.refFile = refFile;
		this.maxRecords = maxRecords;
		this.roundTripCheck = roundTripCheck;
		this.spotAssemblyAlignmentHorizon = spotAssemblyAlignmentHorizon;
		this.spotAssemblyRecordsHorizon = spotAssemblyRecordsHorizon;
		this.maxBlockSize = maxBlockSize;
	}

	public void init() throws IOException {
		samReader = new SAMFileReader(bamFile);
		sequences = new ArrayList<CramReferenceSequence>();
		for (SAMSequenceRecord seq : samReader.getFileHeader()
				.getSequenceDictionary().getSequences()) {
			CramReferenceSequence cramSeq = new CramReferenceSequence(
					seq.getSequenceName(), seq.getSequenceLength());
			sequences.add(cramSeq);
		}
		referenceSequenceFile = ReferenceSequenceFileFactory
				.getReferenceSequenceFile(refFile);
		assembler = new PairedTemplateAssembler(spotAssemblyAlignmentHorizon,
				spotAssemblyRecordsHorizon);

		recordCount = 0;
		unmappedRecordCount = 0;
		baseCount = 0;
		recordsBitLength = 0;
	}

	public void run() throws IOException, CramException {
		for (CramReferenceSequence ref : sequences) {
			if (recordCount >= maxRecords)
				break;
			compressAllRecordsForSequence(ref.getName());
		}
		cramWriter.close();
		if (log.isDebugEnabled())
			cramWriter.dump();

		recordsBitLength = cramWriter.getRecordsBitLength();

		log.info(String.format("Found SAM records: %d\tunmapped: %d",
				recordCount, unmappedRecordCount));
		log.info(String.format("Compressed bases: %d", baseCount));

	}

	private byte[] getReferenceSequenceBases(String seqName) {
		ReferenceSequence sequence = referenceSequenceFile.getSequence(seqName);
		byte[] refBases = referenceSequenceFile.getSubsequenceAt(
				sequence.getName(), 1, sequence.length()).getBases();
		Utils.capitaliseAndCheckBases(refBases);
		return refBases;
	}

	private void compressAllRecordsForSequence(String name) throws IOException,
			CramException {
		assembler.clear();

		cramWriter = new CramWriter(os, provider, sequences, roundTripCheck,
				maxBlockSize);
		cramWriter.init();
		SAMRecordIterator iterator = samReader.query(name, 0, 0, false);
		try {
			if (!iterator.hasNext())
				return;
			byte[] refBases = getReferenceSequenceBases(name);
			cramRecordFactory = new Sam2CramRecordFactory(refBases);
			cramWriter.startSequence(name, refBases);
			while (iterator.hasNext()) {
				if (recordCount >= maxRecords)
					break;

				SAMRecord samRecord = iterator.next();

				addSAMRecord(samRecord);
				recordCount++;

			}
			flushSpotAssembler();
			assembler.clear();
		} finally {
			iterator.close();
		}
	}

	private void flushSpotAssembler() throws IOException, CramException {
		SAMRecord assembledRecord = null;
		while ((assembledRecord = assembler.fetchNextSAMRecord()) != null) {
			writeSAMRecord(assembledRecord, -1);
		}
	}

	private void addSAMRecord(SAMRecord samRecord) throws IOException,
			CramException {
		assembler.addSAMRecord(samRecord);
		SAMRecord assembledRecord = null;
		while ((assembledRecord = assembler.nextSAMRecord()) != null) {
			writeSAMRecord(assembledRecord, assembler.distanceToNextFragment());
		}
	}

	private void writeSAMRecord(SAMRecord record, int distanceToNextFragment)
			throws IOException, CramException {
		CramRecord cramRecord = buildCramRecord(record);
		if (cramRecord.getAlignmentStart() < 1)
			System.out.println("qwerqwer");

		if (!cramRecord.isReadMapped())
			unmappedRecordCount++;

		if (distanceToNextFragment > 0) {
			cramRecord.setLastFragment(false);
			cramRecord.setRecordsToNextFragment(distanceToNextFragment);
		} else
			cramRecord.setLastFragment(true);
		cramWriter.addRecord(cramRecord);
		baseCount += record.getReadLength();
	}

	private CramRecord buildCramRecord(SAMRecord samRecord) {
		return cramRecordFactory.createCramRecord(samRecord);
	}

	private static OutputStream createOutputStream(File outputCramFile,
			boolean wrapInGzip) throws IOException {
		OutputStream os = null;

		if (outputCramFile != null) {
			log.info("Output file: " + outputCramFile.getAbsolutePath());
			FileOutputStream cramFOS = new FileOutputStream(outputCramFile);
			if (wrapInGzip)
				os = new BufferedOutputStream(new GZIPOutputStream(cramFOS));
			else
				os = new BufferedOutputStream(cramFOS);
		} else {
			ByteArrayOutputStream cramBAOS = new ByteArrayOutputStream();
			if (wrapInGzip)
				os = new BufferedOutputStream(new GZIPOutputStream(cramBAOS));
			else
				os = new BufferedOutputStream(cramBAOS);
		}

		return os;
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

		log.info("Input BAM file: " + params.bamFile);
		SAMFileReader
				.setDefaultValidationStringency(ValidationStringency.SILENT);

		OutputStream os = createOutputStream(params.outputFile, params.gzip);

		Bam2Cram bam2Cram_2 = new Bam2Cram(params.bamFile, os,
				params.referenceFasta, params.maxRecords,
				params.roundTripCheck, params.spotAssemblyAlignmentHorizon,
				params.spotAssemblyRecordsHorizon, params.maxBlockSize);
		bam2Cram_2.init();
		long time = System.currentTimeMillis();
		bam2Cram_2.run();
		log.info(String.format("Compression time: %.3f seconds",
				(System.currentTimeMillis() - time) / (float) 1000));
		log.info(String.format("Compression of records: %.4f bits per base.",
				1f * bam2Cram_2.recordsBitLength / bam2Cram_2.baseCount));
		log.info(String.format("Compression, total: %.4f bits per base.", 8f
				* params.outputFile.length() / bam2Cram_2.baseCount));

		os.close();
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

		@Parameter(names = { "--round-trip-check" }, hidden = true)
		boolean roundTripCheck = false;

		@Parameter(names = { "--gzip" }, hidden = true)
		boolean gzip = false;

		@Parameter(names = { "--record-horizon" }, hidden = true)
		int spotAssemblyRecordsHorizon = 10000;

		@Parameter(names = { "--alignment-horizon" }, hidden = true)
		int spotAssemblyAlignmentHorizon = 10000;

		@Parameter(names = { "--max-block-size" }, hidden = true)
		int maxBlockSize = 100000;
	}
}
