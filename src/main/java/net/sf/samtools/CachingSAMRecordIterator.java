package net.sf.samtools;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import net.sf.samtools.util.CloseableIterator;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramReadGroup;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.ReadTag;
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

	Map<String, Integer> seqNameToIndexMap = new TreeMap<String, Integer>();

	public CachingSAMRecordIterator(CloseableIterator<CramRecord> cramRecordIterator, CramHeader cramHeader) {
		this.cramRecordIterator = cramRecordIterator;
		cramReadGroups = cramHeader.getReadGroups();
		samFileHeader = Utils.cramHeader2SamHeader(cramHeader);

		for (SAMSequenceRecord seq : samFileHeader.getSequenceDictionary().getSequences()) {
			seqNameToIndexMap.put(seq.getSequenceName(), seq.getSequenceIndex());
		}

	}

	public CachingSAMRecordIterator(CloseableIterator<CramRecord> cramRecordIterator,
			List<CramReadGroup> cramReadGroups, SAMFileHeader samFileHeader) {
		this.cramRecordIterator = cramRecordIterator;
		this.cramReadGroups = cramReadGroups;
		this.samFileHeader = samFileHeader;
		Map<String, Integer> seqNameToIndexMap = new TreeMap<String, Integer>();
		for (SAMSequenceRecord seq : samFileHeader.getSequenceDictionary().getSequences()) {
			seqNameToIndexMap.put(seq.getSequenceName(), seq.getSequenceIndex());
		}

	}

	private void newSequence(String seqName) {
		currentSeqName = seqName;
		prevAlStart = 0;
		indexes.clear();
	}

	@Override
	public boolean hasNext() {
		boolean hasNext = !assembler.isEmpty() || cramRecordIterator.hasNext();
		return hasNext;
	}

	@Override
	public SAMRecord next() {
		if (cramRecordIterator.hasNext()) {
			SAMRecord samRecord = pollSAMRecordFromAssembler();
			if (samRecord != null) {
				return samRecord;
			}

			while (cramRecordIterator.hasNext()) {
				try {
					readNextCramRecord();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}

				samRecord = pollSAMRecordFromAssembler();
				if (samRecord != null) {
					return samRecord;
				}

			}
		}

		SAMRecord next = takeSAMRecordFromAssembler();
		return next;
	}

	private void fixMateInfo(SAMRecord samRecord) {
		if (!samRecord.getReadPairedFlag()) {
			samRecord.setProperPairFlag(false);
			return;
		}

		if (samRecord.getMateReferenceName() != null
				&& !SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(samRecord.getMateReferenceName()))
			return;
		SAMRecord mate = assembler.getMateRecord();
		if (mate != null)
			Utils.setLooseMateInfo(samRecord, mate, samFileHeader);
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

		SAMRecord samRecord = new SAMRecord(samFileHeader);
		
		if (cramRecord.tags != null && !cramRecord.tags.isEmpty()) {
			for (ReadTag rt : cramRecord.tags) {
				samRecord.setAttribute(rt.getKey(), rt.getValue());
			}
		}

		if (cramReadGroups != null) {
			if (!cramReadGroups.isEmpty()) {
				CramReadGroup cramReadGroup = cramReadGroups.get(cramRecord.getReadGroupID());
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

		samRecord.setMappingQuality(cramRecord.getMappingQuality() & 0xFF);
		if (cramRecord.isReadMapped()) {
			samRecord.setAlignmentStart((int) cramRecord.getAlignmentStart());
			samRecord.setReadBases(cramRecord.getReadBases());
			byte[] scores = cramRecord.getQualityScores();
			injectQualityScores(scores, samRecord);
			prevAlStart = samRecord.getAlignmentStart();
		} else {
			samRecord.setAlignmentStart((int) cramRecord.getAlignmentStart());
			// samRecord.setAlignmentStart((int) prevAlStart);
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

		if (longJump)
			assembler.addSAMRecordNoAssembly(samRecord);
		else
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
