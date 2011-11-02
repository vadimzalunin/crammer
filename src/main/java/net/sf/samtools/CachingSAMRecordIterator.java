package net.sf.samtools;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import net.sf.picard.sam.SamPairUtil;
import net.sf.samtools.util.CloseableIterator;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramReadGroup;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.impl.ReadFeatures2Cigar;
import uk.ac.ebi.ena.sra.cram.spot.PairedTemplateAssembler;

public class CachingSAMRecordIterator implements CloseableIterator<SAMRecord> {
	private static byte defaultQS = '?';

	PairedTemplateAssembler assembler = new PairedTemplateAssembler(10000, 10000);
	CloseableIterator<CramRecord> cramRecordIterator;
	String currentSeqName = null;
	Map<Long, Long> indexes = new TreeMap<Long, Long>();

	long counter = 1;
	long prevAlStart = 0;
	long skipToAlignment = 0;

	List<CramReadGroup> cramReadGroups;
	SAMFileHeader samFileHeader;

	ReadFeatures2Cigar readFeatures2Cigar = new ReadFeatures2Cigar();

	public CachingSAMRecordIterator(CloseableIterator<CramRecord> cramRecordIterator, CramHeader cramHeader) {
		this.cramRecordIterator = cramRecordIterator;
		cramReadGroups = cramHeader.getReadGroups();
		samFileHeader = Utils.cramHeader2SamHeader(cramHeader);
	}

	public CachingSAMRecordIterator(CloseableIterator<CramRecord> cramRecordIterator,
			List<CramReadGroup> cramReadGroups, SAMFileHeader samFileHeader) {
		this.cramRecordIterator = cramRecordIterator;
		this.cramReadGroups = cramReadGroups;
		this.samFileHeader = samFileHeader;
	}

	private void newSequence(String seqName) {
		currentSeqName = seqName;
		prevAlStart = 0;
		indexes.clear();
	}

	@Override
	public boolean hasNext() {
		return !assembler.isEmpty() || cramRecordIterator.hasNext();
	}

	@Override
	public SAMRecord next() {
		if (cramRecordIterator.hasNext()) {
			SAMRecord samRecord = pollSAMRecordFromAssembler();
			if (samRecord != null)
				return samRecord;

			while (cramRecordIterator.hasNext()) {
				try {
					readNextCramRecord();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}

				samRecord = pollSAMRecordFromAssembler();
				if (samRecord != null)
					return samRecord;

			}
		}

		return takeSAMRecordFromAssembler();
	}

	private void fixMateInfo(SAMRecord samRecord) {
		if (!samRecord.getReadPairedFlag()) {
			samRecord.setProperPairFlag(false) ;
			return;
		}

		SAMRecord mate = assembler.getMateRecord();
		if (mate != null) 
			SamPairUtil.setMateInfo(samRecord, mate, samFileHeader);
		else {
			samRecord.setReadPairedFlag(false);
			samRecord.setProperPairFlag(false);
			samRecord.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
			samRecord.setMateNegativeStrandFlag(false);
			samRecord.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
			samRecord.setMateUnmappedFlag(false);
		}
	}

	private SAMRecord takeSAMRecordFromAssembler() {
		SAMRecord samRecord = assembler.fetchNextSAMRecord();
		fixMateInfo(samRecord);
		return samRecord;
	}

	private SAMRecord pollSAMRecordFromAssembler() {
		SAMRecord samRecord = assembler.nextSAMRecord();
		if (samRecord == null)
			return null;

		fixMateInfo(samRecord);
		return samRecord;
	}

	private void readNextCramRecord() throws IOException {
		CramRecord cramRecord = cramRecordIterator.next();
		if (currentSeqName == null || !currentSeqName.equals(cramRecord.getSequenceName()))
			newSequence(cramRecord.getSequenceName());

		if (!cramRecord.isLastFragment()) {
			indexes.put(counter + cramRecord.getRecordsToNextFragment(), counter);
		}

		Long index = indexes.remove(counter);

		SAMRecord samRecord = new SAMRecord(samFileHeader);

		if (cramReadGroups != null) {
			if (!cramReadGroups.isEmpty()) {
				CramReadGroup cramReadGroup = cramReadGroups.get(cramRecord.getReadGroupID());
				String rgId = cramReadGroup.getId();
				if (rgId != null)
					samRecord.setAttribute("RG", rgId);
			}
		}
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

		samRecord.setMappingQuality(cramRecord.getMappingQuality());
		if (cramRecord.isReadMapped()) {
			samRecord.setAlignmentStart((int) cramRecord.getAlignmentStart());
			samRecord.setReadBases(cramRecord.getReadBases());
			byte[] scores = cramRecord.getQualityScores();
			injectQualityScores(scores, samRecord);
			prevAlStart = samRecord.getAlignmentStart();
		} else {
			samRecord.setAlignmentStart((int) prevAlStart);
			samRecord.setReadBases(cramRecord.getReadBases());
			byte[] scores = cramRecord.getQualityScores();
			injectQualityScores(scores, samRecord);
		}
		samRecord
				.setCigar(readFeatures2Cigar.getCigar2(cramRecord.getReadFeatures(), (int) cramRecord.getReadLength()));
		samRecord.setReadUnmappedFlag(!cramRecord.isReadMapped());
		samRecord.setReadNegativeStrandFlag(cramRecord.isNegativeStrand());
		samRecord.setReferenceName(cramRecord.getSequenceName());
		samRecord.setProperPairFlag(cramRecord.isProperPair());
		samRecord.setDuplicateReadFlag(cramRecord.isDuplicate());

		assembler.addSAMRecord(samRecord);

		counter++;
	}

	@Override
	public void remove() {
		throw new RuntimeException("SAM record iterator does not support 'remove' method.");
	}

	@Override
	public void close() {
		indexes.clear();
		assembler.clear();
		cramRecordIterator.close();
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
