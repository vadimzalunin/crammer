package uk.ac.ebi.ena.sra.cram;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMTag;
import net.sf.samtools.util.SequenceUtil;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramReadGroup;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.ReadTag;
import uk.ac.ebi.ena.sra.cram.format.text.CramRecordFormat;
import uk.ac.ebi.ena.sra.cram.impl.ByteArraySequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.impl.CramHeaderIO;
import uk.ac.ebi.ena.sra.cram.impl.ReadFeatures2Cigar;
import uk.ac.ebi.ena.sra.cram.impl.RestoreBases;
import uk.ac.ebi.ena.sra.cram.impl.RestoreQualityScores;
import uk.ac.ebi.ena.sra.cram.impl.SequentialCramReader;
import uk.ac.ebi.ena.sra.cram.spot.PairedTemplateAssembler;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.FileConverter;

public class Cram2Bam {
	private static Logger log = Logger.getLogger(Cram2Bam.class);
	private static byte defaultQS = '?';

	private static InputStream createCramInputStream(File file) throws IOException {
		FileInputStream fis = new FileInputStream(file);

		// gzip magic:
		if (fis.read() == 31 && fis.read() == 139)
			return new GZIPInputStream(new BufferedInputStream(new FileInputStream(file)));

		return new BufferedInputStream(new FileInputStream(file));
	}

	public static void usage(JCommander jc) {
		StringBuilder sb = new StringBuilder();
		sb.append("\n");
		jc.usage(sb);

		System.out.println(sb.toString());
	}

	public static void main(String[] args) throws Exception {
		Params params = new Params();
		JCommander jc = new JCommander(params);
		jc.parse(args);

		if (args.length == 0 || params.help) {
			usage(jc);
			System.exit(1);
		}

		if (params.reference == null) {
			System.out.println("A reference fasta file is required.");
			System.exit(1);
		}

		if (params.cramFile == null) {
			System.out.println("A CRAM file is required.");
			System.exit(1);
		}

		ReferenceSequenceFile referenceSequenceFile = Utils.createIndexedFastaSequenceFile(params.reference);
		InputStream cramIS = null;
		if (params.cramFile == null)
			cramIS = new BufferedInputStream(System.in);
		else
			cramIS = createCramInputStream(params.cramFile);

		// crook nail:
		defaultQS = (byte) params.defaultQS;
		log.info("Using default quality score: " + (char) defaultQS);

		if (params.outputFile == null) {
			if (params.cramFile == null)
				throw new RuntimeException("Output BAM file name is required.");
			params.outputFile = new File(params.cramFile.getAbsolutePath() + ".bam");
		}

		convert(referenceSequenceFile, cramIS, params.outputFile, params.maxRecords, params.printCramRecords,
				params.sequences);

	}

	private static void convert(ReferenceSequenceFile referenceSequenceFile, InputStream cramInputStream,
			File outputBamFile, long maxRecords, boolean printCramRecords, List<String> sequences) throws Exception {

		Utils.isCRAM(cramInputStream);

		DataInputStream cramDIS = new DataInputStream(cramInputStream);

		CramHeader cramHeader = CramHeaderIO.read(Utils.getNextChunk(cramDIS));

		BAMFileWriter writer = new BAMFileWriter(outputBamFile);
		writer.setSortOrder(SortOrder.coordinate, true);
		SAMFileHeader header;

		header = Utils.cramHeader2SamHeader(cramHeader);
		writer.setHeader(header);
		Map<String, Integer> seqNameToIndexMap = new TreeMap<String, Integer>();
		for (SAMSequenceRecord seq : header.getSequenceDictionary().getSequences()) {
			seqNameToIndexMap.put(seq.getSequenceName(), seq.getSequenceIndex());
		}

		long indexOfFirstRecordInTheBlock = 0L;
		CramRecordBlock prevBlock = null;
		long prevAlStart = 0L;

		long time1 = System.currentTimeMillis();
		CramRecordFormat cramRecordFormat = new CramRecordFormat();

		PairedTemplateAssembler assembler = new PairedTemplateAssembler(Integer.MAX_VALUE, Integer.MAX_VALUE);
		String prevSeqName = null;
		Map<Long, Long> indexes = new TreeMap<Long, Long>();

		long counter = 1;
		ByteArraySequenceBaseProvider provider = null;
		byte[] refBases = null;

		NEXT_BLOCK: while (true) {

			DataInputStream nextChunk = Utils.getNextChunk(cramDIS);
			if (nextChunk == null)
				break;

			SequentialCramReader reader = new SequentialCramReader(nextChunk, null, cramHeader);
			CramRecordBlock readBlock = reader.readBlock();
			if (sequences != null && !sequences.contains(readBlock.getSequenceName()))
				continue;
			prevAlStart = readBlock.getFirstRecordPosition();

			if (prevBlock != null) {
				if (!prevBlock.getSequenceName().equals(readBlock.getSequenceName())) {
					prevAlStart = readBlock.getFirstRecordPosition();
				}
			}
			cramRecordFormat.setSequenceID(readBlock.getSequenceName());

			log.info(readBlock);
			SAMSequenceRecord sequence = header.getSequence(readBlock.getSequenceName());
			if (sequence == null) {
				sequence = new SAMSequenceRecord(readBlock.getSequenceName(), readBlock.getSequenceLength());
				header.addSequence(sequence);
			}

			String seqName = readBlock.getSequenceName();
			if (prevSeqName == null)
				prevSeqName = seqName;

			if (!prevSeqName.equals(seqName)) {
				indexes.clear();
				flushAssembler(assembler, writer, header);

				// counter = 1;
			}

			if ("*".equals(seqName)) {
				refBases = new byte[] {};
				provider = new ByteArraySequenceBaseProvider(refBases);
			} else {
				if (provider == null || !seqName.equals(prevSeqName)) {
					refBases = Utils.getReferenceSequenceBases(referenceSequenceFile, seqName);
					provider = new ByteArraySequenceBaseProvider(refBases);
				}
			}
			reader.setReferenceBaseProvider(provider);
			prevSeqName = seqName;

			CramRecord cramRecord = null;
			RestoreBases restoreBases = new RestoreBases();
			restoreBases.setProvider(provider);
			restoreBases.setSequenceName(seqName);
			RestoreQualityScores restoreQualityScores = new RestoreQualityScores();
			ReadFeatures2Cigar readFeatures2Cigar = new ReadFeatures2Cigar();

			ArrayList<CramRecord> records = new ArrayList<CramRecord>((int) readBlock.getRecordCount());
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
			counter -= records.size();

			for (CramRecord record : records) {
				cramRecord = record;

				SAMRecord samRecord = new SAMRecord(header);
				if (record.tags != null && !record.tags.isEmpty()) {
					for (ReadTag rt : record.tags) {
						samRecord.setAttribute(rt.getKey(), rt.getValue());
					}
				}

				if (cramHeader.getReadGroups() != null) {
					if (!cramHeader.getReadGroups().isEmpty()) {
						CramReadGroup cramReadGroup = cramHeader.getReadGroups().get(cramRecord.getReadGroupID());
						String rgId = cramReadGroup.getId();
						if (rgId != null)
							samRecord.setAttribute("RG", rgId);
					}
				}

				boolean longJump = false;
				if (cramRecord.next != null || cramRecord.previous != null) {
					longJump = true;
					CramRecord mate = cramRecord.next == null ? cramRecord.previous : cramRecord.next;
					samRecord.setReadPairedFlag(true);
					samRecord.setMateAlignmentStart((int) mate.getAlignmentStart());
					samRecord.setMateNegativeStrandFlag(mate.isNegativeStrand());
					samRecord.setInferredInsertSize(record.insertSize);

					if (SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(mate.getSequenceName())) {
						samRecord.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
						// samRecord.setMateReferenceName(mate.getSequenceName());
						samRecord.setMateUnmappedFlag(!mate.isReadMapped());
						if (cramRecord.isFirstInPair()) {
							samRecord.setFirstOfPairFlag(true);
							samRecord.setSecondOfPairFlag(false);
						} else {
							samRecord.setFirstOfPairFlag(false);
							samRecord.setSecondOfPairFlag(true);
						}
					} else {
						if (mate.getSequenceName() == null || !seqNameToIndexMap.containsKey(mate.getSequenceName())) {
							samRecord.setReadPairedFlag(false);
							samRecord.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
							samRecord.setMateNegativeStrandFlag(false);
						} else {
							int seqIndex = seqNameToIndexMap.get(mate.getSequenceName());
							samRecord.setMateReferenceIndex(seqIndex);
							samRecord.setMateUnmappedFlag(!mate.isReadMapped());
							if (cramRecord.isFirstInPair()) {
								samRecord.setFirstOfPairFlag(true);
								samRecord.setSecondOfPairFlag(false);
							} else {
								samRecord.setFirstOfPairFlag(false);
								samRecord.setSecondOfPairFlag(true);
							}
						}
					}
					samRecord.setReadName(cramRecord.getReadName());
				} else {
					if (!cramRecord.isLastFragment())
						indexes.put(counter + cramRecord.getRecordsToNextFragment(), counter);

					Long index = indexes.remove(counter);

					if (index != null) {
						samRecord.setReadName(String.valueOf(index.intValue())
						// + ".2"
								);
						samRecord.setFirstOfPairFlag(cramRecord.isFirstInPair());
						samRecord.setSecondOfPairFlag(!cramRecord.isFirstInPair());
						samRecord.setReadPairedFlag(true);
					} else {
						if (cramRecord.isLastFragment()) {
							samRecord.setReadName(String.valueOf(counter));
							samRecord.setReadPairedFlag(false);
							samRecord.setFirstOfPairFlag(false);
							samRecord.setSecondOfPairFlag(false);
						} else {
							samRecord.setReadName(String.valueOf(counter)
							// + ".1"
									);
							samRecord.setFirstOfPairFlag(cramRecord.isFirstInPair());
							samRecord.setSecondOfPairFlag(!cramRecord.isFirstInPair());
							samRecord.setReadPairedFlag(true);
							// samRecord.setMateReferenceName(readBlock.getSequenceName());
						}
					}
				}

				samRecord.setMappingQuality((int) cramRecord.getMappingQuality() & 0xFF);
				if (cramRecord.isReadMapped()) {
					samRecord.setAlignmentStart((int) cramRecord.getAlignmentStart());
					samRecord.setReadBases(restoreBases.restoreReadBases(cramRecord));
					byte[] scores = cramRecord.getQualityScores();
					if (scores == null || scores.length == 0)
						scores = restoreQualityScores.restoreQualityScores(cramRecord);

					injectQualityScores(scores, samRecord);
					prevAlStart = samRecord.getAlignmentStart();
				} else {
					samRecord.setAlignmentStart((int) cramRecord.getAlignmentStart());
					// samRecord.setAlignmentStart((int) prevAlStart);
					samRecord.setReadBases(cramRecord.getReadBases());
					byte[] scores = cramRecord.getQualityScores();
					injectQualityScores(scores, samRecord);
				}
				samRecord.setCigar(readFeatures2Cigar.getCigar2(cramRecord.getReadFeatures(),
						(int) cramRecord.getReadLength()));
				samRecord.setReadUnmappedFlag(!cramRecord.isReadMapped());
				samRecord.setReadNegativeStrandFlag(cramRecord.isNegativeStrand());
				samRecord.setReferenceName(readBlock.getSequenceName());
				samRecord.setProperPairFlag(cramRecord.isProperPair());
				samRecord.setDuplicateReadFlag(cramRecord.isDuplicate());

				if (!samRecord.getReadUnmappedFlag()) 
					Utils.calculateMdAndNmTags(samRecord, refBases);

				if (printCramRecords)
					System.out.println(cramRecordFormat.writeRecord(cramRecord));

				if (longJump)
					assembler.addSAMRecordNoAssembly(samRecord);
				else
					assembler.addSAMRecord(samRecord);

				while ((samRecord = assembler.nextSAMRecord()) != null) {
					SAMRecord mate = assembler.getMateRecord();

					if (mate != null)
						Utils.setLooseMateInfo(samRecord, mate, header);

					writeSAMRecord(samRecord, writer);
				}

				if (counter++ >= maxRecords) {
					indexOfFirstRecordInTheBlock += records.size();
					break NEXT_BLOCK;
				}
			}

			indexOfFirstRecordInTheBlock += records.size();
			prevBlock = readBlock;
		}

		flushAssembler(assembler, writer, header);

		writer.close();
		long time2 = System.currentTimeMillis();
		log.info("Decoded in: " + (time2 - time1) + " millis");
	}

	private static void fixMateInfo(PairedTemplateAssembler assembler, SAMRecord samRecord, SAMFileHeader header) {
		if (!samRecord.getReadPairedFlag()) {
			samRecord.setProperPairFlag(false);
			return;
		}

		if (samRecord.getMateReferenceName() != null
				&& !SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(samRecord.getMateReferenceName()))
			return;
		SAMRecord mate = assembler.getMateRecord();
		if (mate != null)
			Utils.setLooseMateInfo(samRecord, mate, header);
		else {
			if (!samRecord.getReadPairedFlag()) {
				samRecord.setReadPairedFlag(false);
				samRecord.setProperPairFlag(false);
				samRecord.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
				samRecord.setMateNegativeStrandFlag(false);
				samRecord.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
				samRecord.setMateUnmappedFlag(false);
			}
		}
	}

	private static final void flushAssembler(PairedTemplateAssembler assembler, SAMFileWriter writer,
			SAMFileHeader header) {
		SAMRecord samRecord;
		while ((samRecord = assembler.fetchNextSAMRecord()) != null) {
			fixMateInfo(assembler, samRecord, header);
			// SAMRecord mate = assembler.getMateRecord();
			// if (mate != null)
			// SamPairUtil.setMateInfo(samRecord, mate, header);
			writeSAMRecord(samRecord, writer);
		}
	}

	private static final void writeSAMRecord(SAMRecord samRecord, SAMFileWriter writer) {
		try {
			// quick fixes:
			if (SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(samRecord.getReferenceName())) {
				samRecord.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
				samRecord.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
			}

			if (samRecord.getReadUnmappedFlag())
				samRecord.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);
			//

			writer.addAlignment(samRecord);
			// System.out.println(samRecord.format());
		} catch (IllegalArgumentException e) {
			log.error("Offensive SAM record: " + samRecord.format());
			log.error("SAM record al start=" + samRecord.getAlignmentStart());
			throw e;
		}
	}

	@Parameters(commandDescription = "CRAM to BAM conversion. Version 0.65")
	static class Params {
		@Parameter(names = { "--input-cram-file" }, converter = FileConverter.class, description = "The path to the CRAM file to uncompress. Omit if standard input (pipe).")
		File cramFile;

		@Parameter(names = { "--max-sequences" }, description = "Stop after uncomressing this many reference sequences (chromosomes).")
		int maxSequences = Integer.MAX_VALUE;

		@Parameter(names = { "--max-records" }, description = "Stop after uncompressing this many records.")
		long maxRecords = Long.MAX_VALUE;

		@Parameter(names = { "--reference-fasta-file" }, converter = FileConverter.class, description = "Path to the reference fasta file, it must be uncompressed and indexed (use 'samtools faidx' for example).")
		File reference;

		@Parameter(names = { "--output-bam-file" }, converter = FileConverter.class, description = "The path to the output BAM file.")
		File outputFile;

		@Parameter(names = { "-h", "--help" }, description = "Print help and quit")
		boolean help = false;

		@Parameter(names = { "--print-cram-records" }, description = "Print CRAM records while uncompressing.", hidden = true)
		boolean printCramRecords = false;

		@Parameter(names = { "--default-quality-score" }, description = "Use this quality score (decimal representation of ASCII symbol) as a default value when the original quality score was lost due to compression. Minimum is 33.")
		int defaultQS = '?';

		@Parameter()
		List<String> sequences;

	}

	private static final void injectQualityScores(byte[] scores, SAMRecord record) {
		if (scores == null || scores.length == 0) {
			injectNullQualityScores(record);
			return;
		}
		final byte nullQS = -1;
		final byte asciiOffset = 33;
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
				scores[i] = (byte) (defaultQS - asciiOffset);
		}

		record.setBaseQualities(scores);
	}

	private static final void injectNullQualityScores(SAMRecord record) {
		byte[] scores = new byte[record.getReadLength()];
		Arrays.fill(scores, (byte) (defaultQS - 33));
		record.setBaseQualities(scores);
	}
}
