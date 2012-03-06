package uk.ac.ebi.ena.sra.cram;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.zip.GZIPOutputStream;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;

import org.apache.commons.collections.map.MultiValueMap;
import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.bam.Sam2CramRecordFactory;
import uk.ac.ebi.ena.sra.cram.bam.Sam2CramRecordFactory.TREAT_TYPE;
import uk.ac.ebi.ena.sra.cram.format.CramHeaderRecord;
import uk.ac.ebi.ena.sra.cram.format.CramReadGroup;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;
import uk.ac.ebi.ena.sra.cram.format.text.CramRecordFormat;
import uk.ac.ebi.ena.sra.cram.impl.CramWriter;
import uk.ac.ebi.ena.sra.cram.impl.LocalReorderingSAMRecordQueue;
import uk.ac.ebi.ena.sra.cram.impl.ReadAnnotationReader;
import uk.ac.ebi.ena.sra.cram.mask.FastaByteArrayMaskFactory;
import uk.ac.ebi.ena.sra.cram.mask.IntegerListMaskFactory;
import uk.ac.ebi.ena.sra.cram.mask.ReadMaskFactory;
import uk.ac.ebi.ena.sra.cram.mask.RefMaskUtils;
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
	private ReadAnnotationReader readAnnoReader;
	private PrintStream statsPS;
	private PrintStream tramPS;

	private ReferenceSequenceFile referenceSequenceFile;
	private CramRecordFormat cramRecordFormat = new CramRecordFormat();

	private long recordCount = 0;
	private long unmappedRecordCount = 0;
	private long baseCount = 0;
	private long landedRefMaskScores = 0;
	private long landedPiledScores = 0;
	private long landedTotalScores = 0;
	private long beyondHorizon = 0;
	private long extraChromosomePairs = 0;

	private Map<String, Integer> readGroupIdToIndexMap = new TreeMap<String, Integer>();
	private LocalReorderingSAMRecordQueue reoderingQueue = new LocalReorderingSAMRecordQueue(10000);

	private Params params;

	public Bam2Cram(Params params) {
		this.params = params;
	}

	public void init() throws IOException, CramException {
		if (params.bamFile == null) {
			if (params.bamFile == null)
				throw new RuntimeException("Input BAM file name is required.");
		} else
			samReader = new SAMFileReader(params.bamFile);

		samReader.setValidationStringency(ValidationStringency.SILENT);
		sequences = new ArrayList<CramReferenceSequence>();
		referenceSequenceFile = Utils.createIndexedFastaSequenceFile(params.referenceFasta);

		if (!params.exlcudeMappedReads) {
			for (SAMSequenceRecord seq : samReader.getFileHeader().getSequenceDictionary().getSequences()) {
				CramReferenceSequence cramSeq = new CramReferenceSequence(seq.getSequenceName(),
						seq.getSequenceLength());
				sequences.add(cramSeq);
			}
		}

		if (params.inlcudeUnmappedReads)
			sequences.add(new CramReferenceSequence(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME, 0));

		List<CramReadGroup> cramReadGroups = new ArrayList<CramReadGroup>(samReader.getFileHeader().getReadGroups()
				.size() + 1);
		cramReadGroups.add(new CramReadGroup(null));
		for (SAMReadGroupRecord rgr : samReader.getFileHeader().getReadGroups()) {
			readGroupIdToIndexMap.put(rgr.getReadGroupId(), readGroupIdToIndexMap.size() + 1);
			cramReadGroups.add(new CramReadGroup(rgr.getReadGroupId(), rgr.getSample()));
		}

		// copy the whole header:

		assembler = new PairedTemplateAssembler(params.spotAssemblyAlignmentHorizon, params.spotAssemblyRecordsHorizon);

		if (params.readQualityMaskFile != null) {
			log.info("Using read quality mask file: " + params.readQualityMaskFile);
			ReadMaskFactory<String> rqmFactory = params.fastaReadQualityMasking ? new FastaByteArrayMaskFactory()
					: new IntegerListMaskFactory();
			maskReader = new SingleLineMaskReader(new BufferedReader(new FileReader(params.readQualityMaskFile)),
					rqmFactory);
		}

		if (params.readAnnoFile != null) {
			readAnnoReader = new ReadAnnotationReader(new BufferedReader(new FileReader(params.readAnnoFile)));
		}

		recordCount = 0;
		unmappedRecordCount = 0;
		baseCount = 0;

		os = createOutputStream(params.outputCramFile, false);

		statsPS = params.statsOutFile == null ? null : new PrintStream(params.statsOutFile);

		tramPS = params.tramOutFile == null ? null : new PrintStream(params.tramOutFile);

		List<CramHeaderRecord> headerRecords = Utils.getCramHeaderRecords (samReader.getFileHeader()) ;
		cramWriter = new CramWriter(os, provider, sequences, params.roundTripCheck, params.maxBlockSize,
				params.captureUnmappedQualityScore, params.captureSubstitutionQualityScore,
				params.captureMaskedQualityScore, readAnnoReader == null ? null : readAnnoReader.listUniqAnnotations(),
				statsPS, cramReadGroups, params.captureAllQualityScore, headerRecords);
		cramWriter.setAutodump(log.isDebugEnabled());
		cramWriter.init();
	}

	public void close() throws IOException {
		os.close();
		if (statsPS != null)
			statsPS.close();
		if (tramPS != null)
			tramPS.close();

		if (params.outputCramFile != null)
			log.info(String.format("Total compression: %.2f", 8f * params.outputCramFile.length() / baseCount));
	}

	public void run() throws IOException, CramException {
		for (CramReferenceSequence ref : sequences) {
			if (recordCount >= params.maxRecords)
				break;
			if (params.sequences != null && !params.sequences.isEmpty() && !params.sequences.contains(ref.getName()))
				continue;
			compressAllRecordsForSequence2(ref.getName());
		}
		cramWriter.close();

		log.info(String.format("Found SAM records: %d\tunmapped: %d", recordCount, unmappedRecordCount));
		log.info(String.format("Beyond horizon pairs: %d, extra chromosomal pairs: %d", beyondHorizon,
				extraChromosomePairs));
		log.info(String.format("Compressed bases: %d", baseCount));
		log.info(String.format("Landed ref masked qscores: %d", landedRefMaskScores));
		log.info(String.format("Landed piled qscores: %d", landedPiledScores));
		log.info(String.format("Landed total qscores: %d", landedTotalScores));
		log.info(String.format("Quality budget: %.2f%%", baseCount == 0 ? 0 : 100f * landedTotalScores / baseCount));
	}

	private static byte[] refPos2RefSNPs(File file, int refLen, byte onValue) throws IOException {
		if (file == null)
			return null;
		byte[] refSNPs = new byte[refLen];
		BufferedReader r = new BufferedReader(new FileReader(file));
		String line = null;
		while ((line = r.readLine()) != null) {
			int pos = Integer.valueOf(line);
			refSNPs[pos] = onValue;
		}
		r.close();
		return refSNPs;
	}

	private void compressAllRecordsForSequence2(String name) throws IOException, CramException {
		boolean unmapped = SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(name);

		assembler.clear();
		cramRecordFormat.setSequenceID(name);

		SAMRecordIterator iterator;
		if (unmapped)
			iterator = samReader.queryUnmapped();
		else
			iterator = samReader.query(name, 0, 0, false);

		while (params.skipFirstRecords-- > 0 && iterator.hasNext())
			iterator.next();

		long recordsInSequence = 0L;
		try {
			if (!iterator.hasNext())
				return;
			byte[] refBases = null;
			byte[] refSNPs = null;
			RefMaskUtils.RefMask refPileMasks = null;

			if (!unmapped) {
				refBases = Utils.getReferenceSequenceBases(referenceSequenceFile, name);
				refSNPs = refPos2RefSNPs(params.refSnpPosFile, refBases.length, (byte) '+');
				if (params.capturePiledQualityScores) {
					refPileMasks = new RefMaskUtils.RefMask(refBases.length, params.minPiledHits);
				}
			}

			cramRecordFactory = new Sam2CramRecordFactory(refBases, refSNPs, refPileMasks, readGroupIdToIndexMap);
			cramRecordFactory.setCaptureInsertScores(params.captureInsertionQualityScore);
			cramRecordFactory.setCaptureSubtitutionScores(params.captureSubstitutionQualityScore);
			cramRecordFactory.setCaptureUnmappedScores(params.captureUnmappedQualityScore);
			cramRecordFactory.setUncategorisedQualityScoreCutoff(params.qualityCutoff);
			cramRecordFactory.captureAllTags = params.captureAllTags;
			if (params.ignoreTags != null && params.ignoreTags.length() > 0)
				cramRecordFactory.ignoreTags.addAll(Arrays.asList(params.ignoreTags.split(":")));

			if (params.captureAllQualityScore) {
				cramRecordFactory.losslessQS = true;
			} else {
				cramRecordFactory.losslessQS = false;
			}

			if (params.ignoreSoftClips)
				cramRecordFactory.setTreatSoftClipsAs(TREAT_TYPE.IGNORE);
			else
				cramRecordFactory.setTreatSoftClipsAs(TREAT_TYPE.INSERTION);
			boolean enforceAlignmentOrder = false;

			cramWriter.startSequence(name, refBases);

			if (unmapped || !enforceAlignmentOrder) {
				while (iterator.hasNext() && recordCount++ < params.maxRecords
						&& recordsInSequence++ < params.maxRecordsPerSequence) {
					SAMRecord samRecord = iterator.next();
					if ((params.excludeReadsWithFlags & samRecord.getFlags()) != 0 || params.excludeUnmappedPlacedReads
							&& params.excludeUnmappedPlacedReads && samRecord.getReadUnmappedFlag())
						continue;
					compressRecord(samRecord);
				}
			} else {
				/*
				 * this block should take care of the cases when the alignment
				 * order has been changed, for example if soft clips are treated
				 * as match/mismatch.
				 */
				LinkedList<SAMRecord> recordBuffer = new LinkedList<SAMRecord>();
				SAMRecord tempRecord = null;
				int alEnd = 0;
				while (iterator.hasNext() && recordCount++ < params.maxRecords
						&& recordsInSequence++ < params.maxRecordsPerSequence) {
					SAMRecord samRecord = iterator.next();
					if ((params.excludeReadsWithFlags & samRecord.getFlags()) != 0 || params.excludeUnmappedPlacedReads
							&& samRecord.getReadUnmappedFlag())
						continue;

					recordBuffer.add(samRecord);

					if (refPileMasks != null)
						addRefMask(samRecord, refBases, refPileMasks);

					while (!recordBuffer.isEmpty()) {
						tempRecord = recordBuffer.peekFirst();
						if (tempRecord.getReadUnmappedFlag())
							alEnd = tempRecord.getAlignmentStart() + tempRecord.getReadLength();
						else
							alEnd = tempRecord.getAlignmentEnd();

						if (alEnd < samRecord.getAlignmentStart()) {
							compressRecord(recordBuffer.pollFirst());
						} else
							break;
					}
				}
				while (!recordBuffer.isEmpty())
					compressRecord(recordBuffer.pollFirst());
			}

			flushSpotAssembler();
			assembler.clear();

			landedPiledScores += cramRecordFactory.getLandedPiledScores();
			landedRefMaskScores += cramRecordFactory.getLandedRefMaskScores();
			landedTotalScores += cramRecordFactory.getLandedTotalScores();
		} catch (Throwable t) {
			t.printStackTrace();
		} finally {
			iterator.close();
		}
	}

	private void compressRecord(SAMRecord samRecord) throws IOException, CramException {
		if (params.qualityCutoff > 0) {
			byte[] originalScores = new byte[samRecord.getBaseQualities().length];
			System.arraycopy(samRecord.getBaseQualities(), 0, originalScores, 0, originalScores.length);
			for (int i = 0; i < originalScores.length; i++) {
				if (originalScores[i] < params.qualityCutoff)
					samRecord.getBaseQualities()[i] = originalScores[i];
				else
					samRecord.getBaseQualities()[i] = Sam2CramRecordFactory.ignorePositionsWithQualityScore;
			}
		}

		addSAMRecord(samRecord);
	}

	private void addRefMask(SAMRecord record, byte[] refBases, RefMaskUtils.RefMask refMask) {
		int refStartInBlock;
		int readStartInBlock;
		int refStart;
		int readStart;
		byte readBase;
		byte refBase;
		for (AlignmentBlock block : record.getAlignmentBlocks()) {
			refStartInBlock = block.getReferenceStart();
			readStartInBlock = block.getReadStart();
			byte[] readBases = record.getReadBases();
			for (int i = 0; i < block.getLength(); i++) {
				refStart = refStartInBlock + i - 1;
				readStart = readStartInBlock + i - 1;
				readBase = readBases[readStart];
				refBase = refBases[refStart];
				refMask.addReadBase(refStart, readBase, refBase);
			}
		}
	}

	private void flushSpotAssembler() throws IOException, CramException {
		SAMRecord assembledRecord = null;
		while ((assembledRecord = assembler.fetchNextSAMRecord()) != null) {
			writeSAMRecord(assembledRecord, assembler.distanceToNextFragment());
		}
	}

	private static final void rewriteFirstClipIntoInsertion(SAMRecord record) {
		if (record.getAlignmentStart() == record.getUnclippedStart())
			return;

		Cigar cigar = record.getCigar();
		int clipLength = 0;
		LinkedList<CigarElement> newElements = new LinkedList<CigarElement>();
		newElements.addAll(cigar.getCigarElements());

		while (newElements.get(0).getOperator() == CigarOperator.S
				|| newElements.get(0).getOperator() == CigarOperator.H) {
			clipLength += newElements.remove(0).getLength();
		}

		// int elementsToRemove = 0;
		// for (CigarElement ce : cigarElements) {
		// if (ce.getOperator() == CigarOperator.S
		// || ce.getOperator() == CigarOperator.H) {
		// clipLength += ce.getLength();
		// elementsToRemove++;
		// } else
		// break;
		// }

		if (newElements.size() > 1 && newElements.get(1).getOperator() == CigarOperator.I) {
			newElements.removeFirst();

			int newInsertLength = clipLength + newElements.getFirst().getLength();
			newElements.set(0, new CigarElement(newInsertLength, CigarOperator.I));
		} else {
			newElements.addFirst(new CigarElement(clipLength, CigarOperator.M));
		}

		record.setAlignmentStart(record.getUnclippedStart());
		record.setCigar(new Cigar(newElements));
	}

	private void addSAMRecord(SAMRecord samRecord) throws IOException, CramException {
		if (samRecord.getMateReferenceName().equals("="))
			samRecord.setMateReferenceName(samRecord.getReferenceName());
		assembler.addSAMRecord(samRecord);
		SAMRecord assembledRecord = samRecord;
		while ((assembledRecord = assembler.nextSAMRecord()) != null) {
			writeSAMRecord(assembledRecord, assembler.distanceToNextFragment());
		}
	}

	private static final void treatQScores(byte[] scores) {
		for (int i = 0; i < scores.length; i++)
			if (scores[i] == -1)
				scores[i] = '?' - 33;
	}

	private void writeSAMRecord(SAMRecord record, int distanceToNextFragment) throws IOException, CramException {
		CramRecord cramRecord = buildCramRecord(record);
		if (record.getReadGroup() != null) {
			String readGroupId = record.getReadGroup().getReadGroupId();
			Integer readGroupIndex = readGroupIdToIndexMap.get(readGroupId);
			cramRecord.setReadGroupID(readGroupIndex);
		}

		if (readAnnoReader != null)
			cramRecord.setAnnotations(readAnnoReader.nextReadAnnotations());

		if (!cramRecord.isReadMapped())
			unmappedRecordCount++;

		if (record.getReadPairedFlag()) {
			if (distanceToNextFragment > 0) {
				cramRecord.setLastFragment(false);
				cramRecord.setRecordsToNextFragment(distanceToNextFragment);
			} else {
				if (distanceToNextFragment == PairedTemplateAssembler.POINTEE_DISTANCE_NOT_SET)
					cramRecord.setLastFragment(true);
				else {
					cramRecord.detached = true;
					cramRecord.setLastFragment(false);
					if (distanceToNextFragment == PairedTemplateAssembler.DISTANCE_NOT_SET) {
						if (record.getReferenceName().equals(record.getMateReferenceName())) {
							beyondHorizon++;
							// if (Math.abs(record.getInferredInsertSize()) <
							// params.spotAssemblyAlignmentHorizon) {
							// System.out.println("Insert size less then alignment horizon but assemling failed:");
							// assembler.dumpHead();
							// }
						} else
							extraChromosomePairs++;

						cramRecord.setReadName(record.getReadName());
						cramRecord.setSequenceName(record.getReferenceName());
						CramRecord mate = new CramRecord();
						mate.setAlignmentStart(record.getMateAlignmentStart());
						mate.setNegativeStrand(record.getMateNegativeStrandFlag());
						mate.setSequenceName(record.getMateReferenceName());
						mate.setReadName(record.getReadName());
						mate.setReadMapped(!record.getMateUnmappedFlag());

						mate.detached = true;
						if (record.getFirstOfPairFlag()) {
							cramRecord.next = mate;
							mate.previous = cramRecord;
						} else {
							cramRecord.previous = mate;
							mate.next = cramRecord;
						}
					} else
						throw new RuntimeException("Unknown paired distance code: " + distanceToNextFragment);
				}
			}
		} else
			cramRecord.setLastFragment(true);

		cramWriter.addRecord(cramRecord);

		if (tramPS != null)
			tramPS.println(cramRecordFormat.writeRecord(cramRecord));

		if (params.printCramRecords)
			System.out.println(cramRecordFormat.writeRecord(cramRecord));

		// treatQScores (record.getBaseQualities()) ;
		// System.out.println(record.format());
		baseCount += record.getReadLength();

	}

	private CramRecord buildCramRecord(SAMRecord samRecord) {
		return cramRecordFactory.createCramRecord(samRecord);
	}

	private static OutputStream createOutputStream(File outputCramFile, boolean wrapInGzip) throws IOException {
		OutputStream os = null;

		if (outputCramFile != null) {
			FileOutputStream cramFOS = new FileOutputStream(outputCramFile);
			if (wrapInGzip)
				os = new BufferedOutputStream(new GZIPOutputStream(cramFOS));
			else
				os = new BufferedOutputStream(cramFOS);
		} else
			os = new BufferedOutputStream(System.out);

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
			System.out.println("Failed to parse parameteres, detailed message below: ");
			System.out.println(e.getMessage());
			System.out.println();
			System.out.println("See usage: -h");
			System.exit(1);
		}

		if (args.length == 0 || params.help) {
			printUsage(jc);
			System.exit(1);
		}

		if (params.referenceFasta == null) {
			System.out.println("A reference fasta file is required.");
			System.exit(1);
		}

		if (params.bamFile == null) {
			System.out.println("A BAM file is required. ");
			System.exit(1);
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

	@Parameters(commandDescription = "BAM to CRAM converter. Version 0.65")
	static class Params {
		@Parameter(names = { "--input-bam-file" }, converter = FileConverter.class, description = "Path to a BAM file to be converted to CRAM. Omit if standard input (pipe).")
		File bamFile;

		@Parameter(names = { "--reference-fasta-file" }, converter = FileConverter.class, description = "The reference fasta file, uncompressed and indexed (.fai file, use 'samtools faidx'). ")
		File referenceFasta;

		@Parameter(names = { "--output-cram-file" }, converter = FileConverter.class, description = "The path for the output CRAM file. Omit if standard output (pipe).")
		File outputCramFile = null;

		@Parameter(names = { "--max-records" }, description = "Stop after compressing this many records. ")
		long maxRecords = Long.MAX_VALUE;

		@Parameter(names = { "--max-records-per-sequence" }, description = "For each reference sequence (aka chromosome) compress only this many records.")
		long maxRecordsPerSequence = Long.MAX_VALUE;

		@Parameter
		List<String> sequences;

		@Parameter(names = { "-h", "--help" }, description = "Print help and quit")
		boolean help = false;

		@Parameter(names = { "--round-trip-check" }, hidden = true)
		boolean roundTripCheck = false;

		@Parameter(names = { "--record-horizon" }, hidden = true)
		int spotAssemblyRecordsHorizon = 10000;

		@Parameter(names = { "--alignment-horizon" }, hidden = true)
		int spotAssemblyAlignmentHorizon = 10000;

		@Parameter(names = { "--max-block-size" }, hidden = true)
		int maxBlockSize = 100000;

		@Parameter(names = { "--read-quality-mask-file" }, converter = FileConverter.class, description = "Path to the file containing read quality masks.")
		File readQualityMaskFile;

		@Parameter(names = { "--fasta-style-rqm" }, description = "Read quality mask file is in 'fasta' style.")
		boolean fastaReadQualityMasking = false;

		@Parameter(names = { "--capture-unmapped-quality-scores" }, description = "Preserve quality scores for unmapped reads.")
		boolean captureUnmappedQualityScore = false;

		@Parameter(names = { "--capture-substitution-quality-scores" }, description = "Preserve quality scores for substitutions.")
		boolean captureSubstitutionQualityScore = false;

		@Parameter(names = { "--capture-insertion-quality-scores" }, description = "Preserve quality scores for insertions.")
		boolean captureInsertionQualityScore = false;

		@Parameter(names = { "--capture-masked-quality-scores" }, description = "Preserve quality scores for masked bases.")
		boolean captureMaskedQualityScore = false;

		@Parameter(names = { "--capture-all-quality-scores" }, description = "Preserve all quality scores.")
		boolean captureAllQualityScore = false;

		@Parameter(names = { "--print-cram-records" }, description = "Print CRAM records while compressing.")
		boolean printCramRecords = false;

		@Parameter(names = { "--read-anno-file" }, converter = FileConverter.class, description = "Path to the read-level annotations file. ", hidden = true)
		File readAnnoFile;

		@Parameter(names = { "--stats-out-file" }, converter = FileConverter.class, description = "Print detailed statistics to this file.")
		File statsOutFile;

		@Parameter(names = { "--tram-out-file" }, converter = FileConverter.class, description = "Print textual representation of CRAM records to this file.")
		File tramOutFile;

		@Parameter(names = { "--quality-cutoff" }, description = "Preserve quality scores below this value.")
		int qualityCutoff = 0;

		@Parameter(names = { "--ref-snp-pos-file" }, converter = FileConverter.class, description = "Preserve quality scores for positions on the reference listed in this file.")
		File refSnpPosFile;

		@Parameter(names = { "--ignore-soft-clips" }, description = "Treat soft clips as hard clips.")
		boolean ignoreSoftClips = false;

		@Parameter(names = { "--capture-piled-quality-scores" }, description = "Preserve quality score where at least some reads disagree with the reference. See --mini-piled-hits.")
		boolean capturePiledQualityScores = false;

		@Parameter(names = { "--min-piled-hits" }, description = "Preserve quality score where at least this many reads disagree with the reference.")
		int minPiledHits = 2;

		@Parameter(names = { "--include-unmapped-reads" }, description = "Include unmapped reads, namely those that are not assigned to any reference sequence.")
		boolean inlcudeUnmappedReads = false;

		@Parameter(names = { "--exclude-mapped-reads" }, description = "Exclude mapped reads, namely those that are assigned to a reference sequence.")
		boolean exlcudeMappedReads = false;

		@Parameter(names = { "--exclude-unmapped-placed-reads" }, description = "Exclude any unmapped reads, including those with a mapped mate.")
		boolean excludeUnmappedPlacedReads = false;

		@Parameter(names = { "--exclude-reads-with-flags" }, description = "Exclude reads with these bit flags, for example 512 designates the vendor filtered flag.")
		int excludeReadsWithFlags = 0;

		@Parameter(names = { "--capture-all-tags" }, description = "Capture all tags found in the source BAM file.")
		boolean captureAllTags = false;

		@Parameter(names = { "--ignore-tags" }, description = "Do not preserve the tags listed, for example: XA:XO:X0")
		String ignoreTags;

		@Parameter(names = { "--skip-first-records" }, description = "Start compressing records after this many. ", hidden = true)
		int skipFirstRecords = 0;
	}
}
