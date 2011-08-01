package uk.ac.ebi.ena.sra.cram.bam;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import uk.ac.ebi.ena.sra.cram.CramRecordFactory;
import uk.ac.ebi.ena.sra.cram.encoding.BaseChange;
import uk.ac.ebi.ena.sra.cram.format.BaseQualityScore;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.InsertBase;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;

public class Sam2CramRecordFactory implements CramRecordFactory<SAMRecord> {
	private enum TREAT_TYPE {
		IGNORE, ALIGNMENT, INSERTION
	};

	private byte[] refBases;

	/**
	 * Reserved for later use.
	 */
	private TREAT_TYPE treatSoftClipsAs = TREAT_TYPE.INSERTION;

	public final static byte QS_asciiOffset = 33;
	public final static byte unsetQualityScore = 32;
	public final static byte ignorePositionsWithQualityScore = -1;

	public boolean captureUnmappedBases = true;
	public boolean captureUnmappedScores = true;

	public Sam2CramRecordFactory(byte[] refBases) {
		this.refBases = refBases;
	}

	@Override
	public CramRecord createCramRecord(SAMRecord record) {
		CramRecord cramRecord = new CramRecord();
		cramRecord.setAlignmentStart(record.getAlignmentStart());
		cramRecord.setNegativeStrand(record.getReadNegativeStrandFlag());
		cramRecord.setReadMapped(!record.getReadUnmappedFlag()
		// || record.getAlignmentStart() > 0
				);
		cramRecord.setReadLength(record.getReadLength());
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
			List<ReadFeature> features = createVariations(cramRecord, record);
			cramRecord.setReadFeatures(features);
			cramRecord.setPerfectMatch(!record.getReadUnmappedFlag()
					&& features.isEmpty());
		} else
			cramRecord.setPerfectMatch(false);

		if (!cramRecord.isReadMapped()) {
			if (captureUnmappedBases)
				cramRecord.setReadBases(record.getReadBases());
			if (captureUnmappedScores) {
				cramRecord.setQualityScores(record.getBaseQualities());
				for (int i = 0; i < cramRecord.getQualityScores().length; i++)
					cramRecord.getQualityScores()[i] += QS_asciiOffset;
			}
		}

		return cramRecord;
	}

	private boolean isMasked(byte[] scores, int pos) {
		if (scores == null || pos >= scores.length)
			return false;
		return scores[pos] != ignorePositionsWithQualityScore;
	}

	private List<ReadFeature> createVariations(CramRecord cramRecord,
			SAMRecord samRecord) {
		List<ReadFeature> features = new LinkedList<ReadFeature>();
		int zeroBasedPositionInRead = 0;
		int alignmentStartOffset = 0;
		int cigarElementLength = 0;

		List<CigarElement> cigarElements = samRecord.getCigar()
				.getCigarElements();

		byte[] bases = samRecord.getReadBases();
		byte[] qualityScore = samRecord.getBaseQualities();

		for (CigarElement cigarElement : cigarElements) {
			cigarElementLength = cigarElement.getLength();
			CigarOperator operator = cigarElement.getOperator();

			switch (operator) {
			case D:
				features.add(new DeletionVariation(zeroBasedPositionInRead + 1,
						cigarElementLength));
				break;
			case S:
				switch (treatSoftClipsAs) {
				case ALIGNMENT:
					addSubstitutionsAndMaskedBases(cramRecord, features,
							zeroBasedPositionInRead, alignmentStartOffset,
							cigarElementLength, bases, qualityScore);
					break;
				case IGNORE:
					zeroBasedPositionInRead += cigarElementLength;
					break;
				case INSERTION:
					addInsertion(features, zeroBasedPositionInRead,
							cigarElementLength, bases, qualityScore);
					break;
				default:
					throw new IllegalArgumentException(
							"Not sure how to treat soft clips: "
									+ treatSoftClipsAs);
				}
				break;
			case I:
				addInsertion(features, zeroBasedPositionInRead,
						cigarElementLength, bases, qualityScore);
				break;
			case M:
			case X:
			case EQ:
				addSubstitutionsAndMaskedBases(cramRecord, features,
						zeroBasedPositionInRead, alignmentStartOffset,
						cigarElementLength, bases, qualityScore);
				break;
			default:
				throw new IllegalArgumentException(
						"Unsupported cigar operator: "
								+ cigarElement.getOperator());
			}

			if (cigarElement.getOperator().consumesReadBases())
				zeroBasedPositionInRead += cigarElementLength;
			if (cigarElement.getOperator().consumesReferenceBases())
				alignmentStartOffset += cigarElementLength;
		}

		return features;
	}

	private void addInsertion(List<ReadFeature> features,
			int zeroBasedPositionInRead, int cigarElementLength, byte[] bases,
			byte[] scores) {
		byte[] insertedBases = Arrays.copyOfRange(bases,
				zeroBasedPositionInRead, zeroBasedPositionInRead
						+ cigarElementLength);
		for (int i = 0; i < insertedBases.length; i++) {
			InsertBase ib = new InsertBase();
			ib.setPosition(zeroBasedPositionInRead + 1 + i);
			ib.setBase(insertedBases[i]);
			features.add(ib);

			if (isMasked(scores, zeroBasedPositionInRead + i)) {
				byte score = (byte) (QS_asciiOffset + scores[i
						+ zeroBasedPositionInRead]);
				features.add(new BaseQualityScore(zeroBasedPositionInRead + 1
						+ i, score));
			}
		}
	}

	private void addSubstitutionsAndMaskedBases(CramRecord cramRecord,
			List<ReadFeature> features, int fromPosInRead,
			int alignmentStartOffset, int nofReadBases, byte[] bases,
			byte[] qualityScore) {
		int oneBasedPositionInRead;
		for (int i = 0; i < nofReadBases; i++) {
			oneBasedPositionInRead = i + fromPosInRead + 1;

			if (bases[i + fromPosInRead] != refBases[(int) (cramRecord
					.getAlignmentStart() + i + alignmentStartOffset) - 1]) {
				SubstitutionVariation sv = new SubstitutionVariation();
				sv.setPosition(oneBasedPositionInRead);
				sv.setBase(bases[i + fromPosInRead]);
				sv.setRefernceBase(refBases[(int) (cramRecord
						.getAlignmentStart() + i + alignmentStartOffset) - 1]);
				sv.setBaseChange(new BaseChange(sv.getRefernceBase(), sv
						.getBase()));

				features.add(sv);

				if (isMasked(qualityScore, i + fromPosInRead)) {
					byte score = (byte) (QS_asciiOffset + qualityScore[i
							+ fromPosInRead]);
					features.add(new BaseQualityScore(oneBasedPositionInRead,
							score));
				}
			}
		}
	}

	private TREAT_TYPE getTreatSoftClipsAs() {
		return treatSoftClipsAs;
	}

	private void setTreatSoftClipsAs(TREAT_TYPE treatSoftClipsAs) {
		if (treatSoftClipsAs != TREAT_TYPE.INSERTION)
			throw new IllegalArgumentException(
					"Current implemention can treat soft clips as insertions only.");
		this.treatSoftClipsAs = treatSoftClipsAs;
	}

}
