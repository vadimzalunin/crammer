package uk.ac.ebi.ena.sra.cram.bam;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import uk.ac.ebi.ena.sra.cram.CramRecordFactory;
import uk.ac.ebi.ena.sra.cram.encoding.BaseChange;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;

public class Sam2CramRecordFactory implements CramRecordFactory<SAMRecord> {
	public enum TREAT_TYPE {
		IGNORE, ALIGNMENT, INSERTION
	};

	private byte[] refBases;
	private TREAT_TYPE treatSoftClipsAs = TREAT_TYPE.INSERTION;

	public Sam2CramRecordFactory(byte[] refBases) {
		this.refBases = refBases;
	}

	@Override
	public CramRecord createCramRecord(SAMRecord record) {
		CramRecord cramRecord = new CramRecord();
		cramRecord.setAlignmentStart(record.getAlignmentStart());
		cramRecord.setNegativeStrand(record.getReadNegativeStrandFlag());
		cramRecord.setReadMapped(!record.getReadUnmappedFlag() || record.getAlignmentStart() > 0);
		cramRecord.setReadLength(record.getReadLength());
		if (!record.getReadPairedFlag())
			cramRecord.setLastFragment(false);
		else {
			if (record.getFirstOfPairFlag())
				cramRecord.setLastFragment(false);
			else if (record.getSecondOfPairFlag())
				cramRecord.setLastFragment(true);
			else {
//				System.err
//						.println("Could not determine if the read is last fragment in the template: "
//								+ record.toString());
//				throw new IllegalArgumentException(
//						"Could not determine if the read is last fragment in the template: "
//								+ record.toString());
				cramRecord.setLastFragment(true) ;
			}
		}

		List<ReadFeature> features = createVariations(cramRecord, record
				.getCigar().getCigarElements(), record.getReadBases(),
				record.getBaseQualities(), refBases);
		cramRecord.setReadFeatures(features);
		cramRecord.setPerfectMatch(features.isEmpty());

		if (features != null)
			for (ReadFeature f : features) {
				if (f.getPosition() < 1)
					System.out.println(cramRecord);
			}
		
		if (!cramRecord.isReadMapped()) {
			cramRecord.setReadBases(record.getReadBases()) ;
			cramRecord.setQualityScores(record.getBaseQualities()) ;
		}
		return cramRecord;
	}

	private long addjustAlignmentStartToIncludeSoftClips(CramRecord cramRecord,
			List<CigarElement> cigarElements) {
		int posInSeq = 0;
		long readLength = cramRecord.getReadLength();
		for (CigarElement ce : cigarElements) {
			if (ce.getOperator().consumesReferenceBases())
				return cramRecord.getAlignmentStart() + posInSeq;

			if (CigarOperator.SOFT_CLIP == ce.getOperator()) {
				posInSeq -= ce.getLength();
				if (treatSoftClipsAs == TREAT_TYPE.IGNORE)
					readLength -= ce.getLength();
			}
		}
		cramRecord.setReadLength(readLength);
		throw new IllegalArgumentException(
				"Reference bases are not consumed in this cigar: "
						+ new Cigar(cigarElements));
	}

	private List<ReadFeature> createVariations(CramRecord cramRecord,
			List<CigarElement> cigarElements, byte[] bases,
			byte[] qualityScore, byte[] refBases) {
		List<ReadFeature> features = new LinkedList<ReadFeature>();
		int posInRead = 0;
		int posInSeq = 0;
		int ceLen = 0;

		for (CigarElement ce : cigarElements) {
			ceLen = ce.getLength();
			CigarOperator operator = ce.getOperator();

			switch (operator) {
			case D:
				int pos = addjustPositionIfNeeded(posInRead,
						cramRecord.getReadLength(),
						cramRecord.isNegativeStrand());
				features.add(new DeletionVariation(pos, ceLen));
				posInSeq += ceLen;
				break;
			case S:
				switch (treatSoftClipsAs) {
				case ALIGNMENT:
					for (int i = 0; i < ceLen; i++) {
						if (bases[posInRead] != refBases[(int) (cramRecord
								.getAlignmentStart() + posInSeq) - 1]) {
							SubstitutionVariation sv = new SubstitutionVariation();
							pos = addjustPositionIfNeeded(posInRead,
									cramRecord.getReadLength(),
									cramRecord.isNegativeStrand());
							sv.setPosition(pos);
							sv.setBase(bases[posInRead]);
							sv.setRefernceBase(refBases[(int) (cramRecord
									.getAlignmentStart() + posInSeq) - 1]);
							sv.setBaseChange(new BaseChange(sv
									.getRefernceBase(), sv.getBase()));
							sv.setQualityScore(qualityScore[posInRead]);
							features.add(sv);
						}
						posInRead++;
						posInSeq++;
					}
					break;
				case IGNORE:
					posInRead += ceLen;
					break;
				case INSERTION:
					byte[] insertedBases = Arrays.copyOfRange(bases, posInRead,
							posInRead + ceLen);
					pos = addjustPositionIfNeeded(posInRead,
							cramRecord.getReadLength(),
							cramRecord.isNegativeStrand());
					features.add(new InsertionVariation(pos, insertedBases));
					posInRead += ceLen;
					break;
				default:
					throw new IllegalArgumentException(
							"Not sure how to treat soft clips: "
									+ treatSoftClipsAs);
				}
				break;
			case I:
				byte[] insertedBases = Arrays.copyOfRange(bases, posInRead,
						posInRead + ceLen);
				pos = addjustPositionIfNeeded(posInRead,
						cramRecord.getReadLength(),
						cramRecord.isNegativeStrand());
				features.add(new InsertionVariation(pos, insertedBases));
				posInRead += ceLen;
				break;
			case M:
			case X:
			case EQ:
				for (int i = 0; i < ceLen; i++) {
					if (bases[posInRead] != refBases[(int) (cramRecord
							.getAlignmentStart() + posInSeq) - 1]) {
						SubstitutionVariation sv = new SubstitutionVariation();
						pos = addjustPositionIfNeeded(posInRead,
								cramRecord.getReadLength(),
								cramRecord.isNegativeStrand());
						sv.setPosition(pos);
						sv.setBase(bases[posInRead]);
						sv.setRefernceBase(refBases[(int) (cramRecord
								.getAlignmentStart() + posInSeq) - 1]);
						sv.setBaseChange(new BaseChange(sv.getRefernceBase(),
								sv.getBase()));
						sv.setQualityScore(qualityScore[posInRead]);
						features.add(sv);
					}
					posInRead++;
					posInSeq++;
				}
				break;
			default:
				throw new IllegalArgumentException(
						"Unsupported cigar operator: " + ce.getOperator());
			}
		}

		return features;
	}

	@Deprecated
	private final int addjustPositionIfNeeded(int pos, long legnth,
			boolean negativeStrand) {
		return pos + 1;
	}

	public TREAT_TYPE getTreatSoftClipsAs() {
		return treatSoftClipsAs;
	}

	public void setTreatSoftClipsAs(TREAT_TYPE treatSoftClipsAs) {
		if (treatSoftClipsAs != TREAT_TYPE.INSERTION)
			throw new IllegalArgumentException(
					"Current implemention can treat soft clips as insertions only.");
		this.treatSoftClipsAs = treatSoftClipsAs;
	}
}
