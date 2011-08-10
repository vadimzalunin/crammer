package uk.ac.ebi.ena.sra.cram;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.bam.SAMUtils;
import uk.ac.ebi.ena.sra.cram.bam.Sam2CramRecordFactory;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;
import uk.ac.ebi.ena.sra.cram.format.text.CramRecordFormat;
import uk.ac.ebi.ena.sra.cram.impl.CramWriter;
//import uk.ac.ebi.ena.sra.cram.impl.ReadAnnotationReader;
import uk.ac.ebi.ena.sra.cram.mask.FastaByteArrayMaskFactory;
import uk.ac.ebi.ena.sra.cram.mask.IntegerListMaskFactory;
import uk.ac.ebi.ena.sra.cram.mask.PositionMask;
import uk.ac.ebi.ena.sra.cram.mask.ReadMaskFactory;
import uk.ac.ebi.ena.sra.cram.mask.SingleLineMaskReader;
import uk.ac.ebi.ena.sra.cram.spot.PairedTemplateAssembler;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.FileConverter;

public class Bam2Cram {

	private static Logger log = Logger.getLogger(Bam2Cram.class);

	private SAMFileReader samReader;
	private CramWriter cramWriter;
	private PairedTemplateAssembler assembler;
	private Sam2CramRecordFactory cramRecordFactory;
	private List<CramReferenceSequence> sequences;
	private OutputStream os;
	private SequenceBaseProvider provider;
	private SingleLineMaskReader maskReader;
	// private ReadAnnotationReader readAnnoReader;

	private ReferenceSequenceFile referenceSequenceFile;
	private CramRecordFormat cramRecordFormat = new CramRecordFormat();

	private long recordCount;
	private long unmappedRecordCount;
	private long baseCount;

	private Params params;

	public Bam2Cram(Params params) {
		this.params = params;
	}

	public void init() throws IOException, CramException {
		log.info("Input BAM file: " + params.bamFile.getAbsolutePath());

		samReader = new SAMFileReader(params.bamFile);
		samReader.setValidationStringency(ValidationStringency.SILENT);
		sequences = new ArrayList<CramReferenceSequence>();
		for (SAMSequenceRecord seq : samReader.getFileHeader()
				.getSequenceDictionary().getSequences()) {
			if (params.sequences != null && !params.sequences.isEmpty()
					&& !params.sequences.contains(seq.getSequenceName()))
				continue;
			CramReferenceSequence cramSeq = new CramReferenceSequence(
					seq.getSequenceName(), seq.getSequenceLength());
			sequences.add(cramSeq);
		}
		referenceSequenceFile = SAMUtils
				.createIndexedFastaSequenceFile(params.referenceFasta);
		assembler = new PairedTemplateAssembler(
				params.spotAssemblyAlignmentHorizon,
				params.spotAssemblyRecordsHorizon);

		if (params.readQualityMaskFile != null) {
			log.info("Using read quality mask file: "
					+ params.readQualityMaskFile);
			ReadMaskFactory<String> rqmFactory = params.fastaReadQualityMasking ? new FastaByteArrayMaskFactory()
					: new IntegerListMaskFactory();
			maskReader = new SingleLineMaskReader(new BufferedReader(
					new FileReader(params.readQualityMaskFile)), rqmFactory);
		}

		// if (params.readAnnoFile != null) {
		// readAnnoReader = new ReadAnnotationReader(new BufferedReader(
		// new FileReader(params.readAnnoFile)));
		// }

		recordCount = 0;
		unmappedRecordCount = 0;
		baseCount = 0;

		if (params.outputCramFile != null)
			log.info("Output CRAM file: "
					+ params.outputCramFile.getAbsolutePath());
		else
			log.info("No output CRAM file specified, discarding CRAM output.");

		os = createOutputStream(params.outputCramFile,
				params.gzipOutputCramFile);
		cramWriter = new CramWriter(os, provider, sequences,
				params.roundTripCheck, params.maxBlockSize,
				params.captureUnmappedQualityScore,
				params.captureSubstitutionQualityScore,
				params.captureMaskedQualityScore, null);
		cramWriter.setAutodump(log.isInfoEnabled());
		cramWriter.init();
	}

	public void close() throws IOException {
		os.close();
	}

	public void run() throws IOException, CramException {
		for (CramReferenceSequence ref : sequences) {
			if (recordCount >= params.maxRecords)
				break;
			compressAllRecordsForSequence(ref.getName());
		}
		cramWriter.close();

		log.info(String.format("Found SAM records: %d\tunmapped: %d",
				recordCount, unmappedRecordCount));
		log.info(String.format("Compressed bases: %d", baseCount));

	}

	private byte[] getReferenceSequenceBases(String seqName) {
		ReferenceSequence sequence = referenceSequenceFile.getSequence(seqName);
		byte[] refBases = referenceSequenceFile.getSubsequenceAt(
				sequence.getName(), 1, sequence.length()).getBases();
		Utils.capitaliseAndCheckBases(refBases, true);
		return refBases;
	}

	private void compressAllRecordsForSequence(String name) throws IOException,
			CramException {
		assembler.clear();

		SAMRecordIterator iterator = samReader.query(name, 0, 0, false);

		long recordsInSequence = 0L;
		try {
			if (!iterator.hasNext())
				return;
			byte[] refBases = getReferenceSequenceBases(name);
			cramRecordFactory = new Sam2CramRecordFactory(refBases);
			cramWriter.startSequence(name, refBases);
			while (iterator.hasNext()) {
				if (recordCount >= params.maxRecords)
					break;
				if (recordsInSequence >= params.maxRecordsPerSequence)
					break;

				SAMRecord samRecord = iterator.next();
				PositionMask mask = null;
				if (maskReader != null) {
					mask = maskReader.readNextMask();
					byte[] scores = samRecord.getBaseQualities();
					if (mask == null) {
						Arrays.fill(
								samRecord.getBaseQualities(),
								Sam2CramRecordFactory.ignorePositionsWithQualityScore);
					} else if (!mask.isEmpty()) {
						for (int i = 0; i < scores.length; i++)
							if (!mask.isMasked(i + 1))
								scores[i] = Sam2CramRecordFactory.ignorePositionsWithQualityScore;
					} else {
						Arrays.fill(
								samRecord.getBaseQualities(),
								Sam2CramRecordFactory.ignorePositionsWithQualityScore);
					}
				} else
					Arrays.fill(
							samRecord.getBaseQualities(),
							Sam2CramRecordFactory.ignorePositionsWithQualityScore);

				addSAMRecord(samRecord);
				recordCount++;
				recordsInSequence++;

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

		// if (readAnnoReader != null)
		// cramRecord.setAnnotations(readAnnoReader.nextReadAnnotations());

		if (!cramRecord.isReadMapped())
			unmappedRecordCount++;

		if (distanceToNextFragment > 0) {
			cramRecord.setLastFragment(false);
			cramRecord.setRecordsToNextFragment(distanceToNextFragment);
		} else
			cramRecord.setLastFragment(true);
		cramWriter.addRecord(cramRecord);
		if (params.printCramRecords)
			System.out.println(cramRecordFormat.writeRecord(cramRecord));
		baseCount += record.getReadLength();
	}

	private CramRecord buildCramRecord(SAMRecord samRecord) {
		return cramRecordFactory.createCramRecord(samRecord);
	}

	private static OutputStream createOutputStream(File outputCramFile,
			boolean wrapInGzip) throws IOException {
		OutputStream os = null;

		if (outputCramFile != null) {
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

	private static void printUsage(JCommander jc) {
		StringBuilder sb = new StringBuilder();
		sb.append("\n");
		jc.usage(sb);

		System.out.println(sb.toString());
	}

	public static void main(String[] args) throws Exception {
		Params params = new Params();
		JCommander jc = new JCommander(params);
		try {
			jc.parse(args);
		} catch (Exception e) {
			System.out.println(e.getMessage());
			printUsage(jc);
			return;
		}

		if (args.length == 0 || params.help) {
			printUsage(jc);
			return;
		}

		Bam2Cram b2c = new Bam2Cram(params);
		b2c.init();
		b2c.run();
		b2c.close();

		// long time = System.currentTimeMillis();
		// log.info(String.format("Compression time: %.3f seconds",
		// (System.currentTimeMillis() - time) / (float) 1000));
		// if (params.outputFile != null)
		// log.info(String.format("Compression, total: %.4f bits per base.",
		// 8f * params.outputFile.length() / bam2Cram_2.baseCount));
		//
		// os.close();
	}

	@Parameters(commandDescription = "BAM to CRAM converter.")
	static class Params {
		@Parameter(names = { "--input-bam-file" }, converter = FileConverter.class, required = true)
		File bamFile;

		@Parameter(names = { "--reference-fasta-file" }, converter = FileConverter.class, required = true)
		File referenceFasta;

		@Parameter(names = { "--output-cram-file" }, converter = FileConverter.class)
		File outputCramFile = null;

		@Parameter(names = { "--max-records" })
		long maxRecords = Long.MAX_VALUE;

		@Parameter(names = { "--max-records-per-sequence" })
		long maxRecordsPerSequence = Long.MAX_VALUE;

		@Parameter
		List<String> sequences;

		@Parameter(names = { "-h", "--help" }, description = "Print help and quit")
		boolean help = false;

		@Parameter(names = { "--round-trip-check" }, hidden = true)
		boolean roundTripCheck = false;

		@Parameter(names = { "--gzip" })
		boolean gzipOutputCramFile = false;

		@Parameter(names = { "--record-horizon" }, hidden = true)
		int spotAssemblyRecordsHorizon = 10000;

		@Parameter(names = { "--alignment-horizon" }, hidden = true)
		int spotAssemblyAlignmentHorizon = 10000;

		@Parameter(names = { "--max-block-size" }, hidden = true)
		int maxBlockSize = 100000;

		@Parameter(names = { "--read-quality-mask-file" }, converter = FileConverter.class)
		File readQualityMaskFile;

		@Parameter(names = { "--fasta-style-rqm" })
		boolean fastaReadQualityMasking = false;

		@Parameter(names = { "--capture-unmapped-quality-scores" }, description = "Unused", hidden = true)
		boolean captureUnmappedQualityScore = false;

		@Parameter(names = { "--capture-substitution-quality-scores" }, description = "Unused", hidden = true)
		boolean captureSubstitutionQualityScore = false;

		@Parameter(names = { "--capture-masked-quality-scores" }, description = "Unused", hidden = true)
		boolean captureMaskedQualityScore = false;

		@Parameter(names = { "--print-cram-records" })
		boolean printCramRecords = false;

		@Parameter(names = { "--read-anno-file" }, converter = FileConverter.class, hidden = false)
		File readAnnoFile;
	}
}
