package uk.ac.ebi.ena.sra.cram.bam;

import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.CramRecordFactory;
import uk.ac.ebi.ena.sra.cram.encoding.BaseChange;
import uk.ac.ebi.ena.sra.cram.format.BaseQualityScore;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.InsertBase;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;
import uk.ac.ebi.ena.sra.cram.mask.RefMaskUtils;

public class Sam2CramRecordFactory implements CramRecordFactory<SAMRecord> {
	public enum TREAT_TYPE {
		IGNORE, ALIGNMENT, INSERTION
	};

	/**
	 * Reserved for later use.
	 */
	private TREAT_TYPE treatSoftClipsAs = TREAT_TYPE.IGNORE;

	public final static byte QS_asciiOffset = 33;
	public final static byte unsetQualityScore = 32;
	public final static byte ignorePositionsWithQualityScore = -1;

	private byte[] refBases;
	private byte[] refSNPs;
	private RefMaskUtils.RefMask refPile;

	public boolean captureUnmappedBases = true;
	public boolean captureUnmappedScores = false;

	private ByteBuffer insertionBuf = ByteBuffer.allocate(1024);

	private static Logger log = Logger.getLogger(Sam2CramRecordFactory.class);

	private Map<String, Integer> readGroupMap;

	private long landedRefMaskScores = 0;
	private long landedPiledScores = 0;
	private long landedTotalScores = 0;

	private boolean captureInsertScores = false;
	private boolean captureSubtitutionScores = false;
	private boolean captureFlankingDeletionScores = false;
	private int uncategorisedQualityScoreCutoff = 0;
	
	public boolean losslessQS = false;

	public Sam2CramRecordFactory(byte[] refBases) {
		this(refBases, null, null, null);
	}

	public Sam2CramRecordFactory(byte[] refBases, byte[] refSNPs, RefMaskUtils.RefMask refPile,
			Map<String, Integer> readGroupMap) {
		this.refPile = refPile;
//		if (refBases == null)
//			throw new NullPointerException("Reference bases array is null.");
		this.refBases = refBases;
		this.refSNPs = refSNPs;
		this.readGroupMap = readGroupMap;
	}

	public Sam2CramRecordFactory() {
	}

	@Override
	public CramRecord createCramRecord(SAMRecord record) {
		CramRecord cramRecord = new CramRecord();
		cramRecord.setAlignmentStart(record.getAlignmentStart());
		cramRecord.setNegativeStrand(record.getReadNegativeStrandFlag());
		cramRecord.setReadMapped(!record.getReadUnmappedFlag());
		cramRecord.setReadLength(record.getReadLength());
		cramRecord.setFirstInPair(record.getReadPairedFlag() && record.getFirstOfPairFlag());

		cramRecord.setProperPair(record.getReadPairedFlag() && record.getProperPairFlag());
		cramRecord.setMappingQuality((byte) record.getMappingQuality());
		cramRecord.setDuplicate(record.getDuplicateReadFlag());

		if (readGroupMap != null) {
			SAMReadGroupRecord readGroup = record.getReadGroup();
			if (readGroup != null) {
				Integer rgIndex = readGroupMap.get(readGroup.getId());
				if (rgIndex == null)
					throw new RuntimeException("Read group index not found: " + readGroup.getId());
				cramRecord.setReadGroupID(rgIndex);
			}
		}
		if (!record.getReadPairedFlag())
			cramRecord.setLastFragment(false);
		else {
			if (record.getFirstOfPairFlag())
				cramRecord.setLastFragment(false);
			else if (record.getSecondOfPairFlag())
				cramRecord.setLastFragment(true);
			else
				cramRecord.setLastFragment(true);
		}

		if (cramRecord.isReadMapped()) {
			List<ReadFeature> features = checkedCreateVariations(cramRecord, record);
			cramRecord.setReadFeatures(features);
			cramRecord.setPerfectMatch(!record.getReadUnmappedFlag() && features.isEmpty());
		} else {
			cramRecord.setPerfectMatch(false);
			// if (captureUnmappedBases)
			// cramRecord.setReadBases(record.getReadBases());
			// if (captureUnmappedScores) {
			// cramRecord.setQualityScores(record.getBaseQualities());
			// for (int i = 0; i < cramRecord.getQualityScores().length; i++)
			// cramRecord.getQualityScores()[i] += QS_asciiOffset;
			// landedTotalScores += cramRecord.getReadLength();
			// }
		}

		cramRecord.setReadBases(record.getReadBases());
		cramRecord.setQualityScores(record.getBaseQualities());
		for (int i = 0; i < cramRecord.getQualityScores().length; i++)
			cramRecord.getQualityScores()[i] += QS_asciiOffset;
		landedTotalScores += cramRecord.getReadLength();

		return cramRecord;
	}

	/**
	 * A wrapper method to provide better diagnostics for
	 * ArrayIndexOutOfBoundsException.
	 * 
	 * @param cramRecord
	 * @param samRecord
	 * @return
	 */
	private List<ReadFeature> checkedCreateVariations(CramRecord cramRecord, SAMRecord samRecord) {
		try {
			return createVariations(cramRecord, samRecord);
		} catch (ArrayIndexOutOfBoundsException e) {
			log.error("Reference bases array length=" + refBases.length);
			log.error("Offensive CRAM record: " + cramRecord.toString());
			log.error("Offensive SAM record: " + samRecord.getSAMString());
			throw e;
		}
	}

	private List<ReadFeature> createVariations(CramRecord cramRecord, SAMRecord samRecord) {
		List<ReadFeature> features = new LinkedList<ReadFeature>();
		int zeroBasedPositionInRead = 0;
		int alignmentStartOffset = 0;
		int cigarElementLength = 0;

		List<CigarElement> cigarElements = samRecord.getCigar().getCigarElements();

		byte[] bases = samRecord.getReadBases();
		byte[] qualityScore = samRecord.getBaseQualities();

		for (CigarElement cigarElement : cigarElements) {
			cigarElementLength = cigarElement.getLength();
			CigarOperator operator = cigarElement.getOperator();

			switch (operator) {
			case D:
			case N:
				features.add(new DeletionVariation(zeroBasedPositionInRead + 1, cigarElementLength));
				break;
			case S:
				switch (treatSoftClipsAs) {
				// case ALIGNMENT:
				// addSubstitutionsAndMaskedBases(cramRecord, features,
				// zeroBasedPositionInRead, alignmentStartOffset,
				// cigarElementLength, bases, qualityScore);
				// break;
				case IGNORE:
					// zeroBasedPositionInRead += cigarElementLength;
					break;
				case INSERTION:
					addInsertion(features, zeroBasedPositionInRead, cigarElementLength, bases, qualityScore);
					break;
				default:
					throw new IllegalArgumentException("Not sure how to treat soft clips: " + treatSoftClipsAs);
				}
				break;
			case I:
				addInsertion(features, zeroBasedPositionInRead, cigarElementLength, bases, qualityScore);
				break;
			case M:
			case X:
			case EQ:
				addSubstitutionsAndMaskedBases(cramRecord, features, zeroBasedPositionInRead, alignmentStartOffset,
						cigarElementLength, bases, qualityScore);
				break;
			default:
				throw new IllegalArgumentException("Unsupported cigar operator: " + cigarElement.getOperator());
			}

			if (cigarElement.getOperator().consumesReadBases())
				zeroBasedPositionInRead += cigarElementLength;
			if (cigarElement.getOperator().consumesReferenceBases())
				alignmentStartOffset += cigarElementLength;
		}

		return features;
	}

	private void addInsertion(List<ReadFeature> features, int zeroBasedPositionInRead, int cigarElementLength,
			byte[] bases, byte[] scores) {
		byte[] insertedBases = Arrays.copyOfRange(bases, zeroBasedPositionInRead, zeroBasedPositionInRead
				+ cigarElementLength);

		for (int i = 0; i < insertedBases.length; i++) {
			// single base insertion:
			InsertBase ib = new InsertBase();
			ib.setPosition(zeroBasedPositionInRead + 1 + i);
			ib.setBase(insertedBases[i]);
			features.add(ib);
			if (losslessQS) continue ;
			boolean qualityMasked = (scores[i] < uncategorisedQualityScoreCutoff);
			if (captureInsertScores || qualityMasked) {
				byte score = (byte) (QS_asciiOffset + scores[zeroBasedPositionInRead + i]);
				// if (score >= QS_asciiOffset) {
				features.add(new BaseQualityScore(zeroBasedPositionInRead + 1 + i, score));
				landedTotalScores++;
				// }
			}
		}
	}

	private void addSubstitutionsAndMaskedBases(CramRecord cramRecord, List<ReadFeature> features, int fromPosInRead,
			int alignmentStartOffset, int nofReadBases, byte[] bases, byte[] qualityScore) {
		int oneBasedPositionInRead;

		int i = 0;
		boolean qualityAdded = false;
		boolean qualityMasked = false;
		for (i = 0; i < nofReadBases; i++) {
			oneBasedPositionInRead = i + fromPosInRead + 1;
			int refCoord = (int) (cramRecord.getAlignmentStart() + i + alignmentStartOffset) - 1;
			qualityAdded = false;

			if (bases[i + fromPosInRead] != refBases[refCoord]) {
				SubstitutionVariation sv = new SubstitutionVariation();
				sv.setPosition(oneBasedPositionInRead);
				sv.setBase(bases[i + fromPosInRead]);
				sv.setRefernceBase(refBases[refCoord]);
				sv.setBaseChange(new BaseChange(sv.getRefernceBase(), sv.getBase()));

				features.add(sv);
				
				if (losslessQS) continue ;

				if (captureSubtitutionScores) {
					byte score = (byte) (QS_asciiOffset + qualityScore[i + fromPosInRead]);
					features.add(new BaseQualityScore(oneBasedPositionInRead, score));
					qualityAdded = true;
				}
			}
			if (!qualityAdded && refSNPs != null) {
				byte snpOrNot = refSNPs[refCoord];
				if (snpOrNot != 0) {
					byte score = (byte) (QS_asciiOffset + qualityScore[i + fromPosInRead]);
					features.add(new BaseQualityScore(oneBasedPositionInRead, score));
					qualityAdded = true;
					landedRefMaskScores++;
				}
			}

			if (!qualityAdded && refPile != null) {
				if (refPile.shouldStore(refCoord, refBases[refCoord])) {
					byte score = (byte) (QS_asciiOffset + qualityScore[i + fromPosInRead]);
					features.add(new BaseQualityScore(oneBasedPositionInRead, score));
					qualityAdded = true;
					landedPiledScores++;
				}
			}

			qualityMasked = (qualityScore[i + fromPosInRead] < uncategorisedQualityScoreCutoff);
			if (!qualityAdded && qualityMasked) {
				byte score = (byte) (QS_asciiOffset + qualityScore[i + fromPosInRead]);
				features.add(new BaseQualityScore(oneBasedPositionInRead, score));
				qualityAdded = true;
			}

			if (qualityAdded)
				landedTotalScores++;
		}
	}

	public TREAT_TYPE getTreatSoftClipsAs() {
		return treatSoftClipsAs;
	}

	public void setTreatSoftClipsAs(TREAT_TYPE treatSoftClipsAs) {
		if (treatSoftClipsAs == TREAT_TYPE.ALIGNMENT)
			throw new IllegalArgumentException("Current implemention cannot treat soft clips as alignmen, sorry.");
		this.treatSoftClipsAs = treatSoftClipsAs;
	}

	public boolean isCaptureInsertScores() {
		return captureInsertScores;
	}

	public void setCaptureInsertScores(boolean captureInsertScores) {
		this.captureInsertScores = captureInsertScores;
	}

	public boolean isCaptureSubtitutionScores() {
		return captureSubtitutionScores;
	}

	public void setCaptureSubtitutionScores(boolean captureSubtitutionScores) {
		this.captureSubtitutionScores = captureSubtitutionScores;
	}

	public int getUncategorisedQualityScoreCutoff() {
		return uncategorisedQualityScoreCutoff;
	}

	public void setUncategorisedQualityScoreCutoff(int uncategorisedQualityScoreCutoff) {
		this.uncategorisedQualityScoreCutoff = uncategorisedQualityScoreCutoff;
	}

	public long getLandedRefMaskScores() {
		return landedRefMaskScores;
	}

	public long getLandedPiledScores() {
		return landedPiledScores;
	}

	public long getLandedTotalScores() {
		return landedTotalScores;
	}

	public boolean isCaptureUnmappedBases() {
		return captureUnmappedBases;
	}

	public void setCaptureUnmappedBases(boolean captureUnmappedBases) {
		this.captureUnmappedBases = captureUnmappedBases;
	}

	public boolean isCaptureUnmappedScores() {
		return captureUnmappedScores;
	}

	public void setCaptureUnmappedScores(boolean captureUnmappedScores) {
		this.captureUnmappedScores = captureUnmappedScores;
	}

	public byte[] getRefBases() {
		return refBases;
	}

	public void setRefBases(byte[] refBases) {
		this.refBases = refBases;
	}

	public byte[] getRefSNPs() {
		return refSNPs;
	}

	public void setRefSNPs(byte[] refSNPs) {
		this.refSNPs = refSNPs;
	}

	public RefMaskUtils.RefMask getRefPile() {
		return refPile;
	}

	public Map<String, Integer> getReadGroupMap() {
		return readGroupMap;
	}

	public void setRefPile(RefMaskUtils.RefMask refPile) {
		this.refPile = refPile;
	}

}
