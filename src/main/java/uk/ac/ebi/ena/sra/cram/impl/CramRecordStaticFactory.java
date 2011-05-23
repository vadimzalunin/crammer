package uk.ac.ebi.ena.sra.cram.impl;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;

public class CramRecordStaticFactory {

	private static void addVariations(CramRecord cramRecord,
			List<CigarElement> cigarElements, byte[] bases,
			byte[] qualityScore, byte[] refBases) {
		long readLength = cramRecord.getReadLength();
		int posInRead = 0;
		int posInSeq = 0;
		boolean perfectMatch = true;
		int ceLen = 0;
		for (CigarElement ce : cigarElements) {

			ceLen = ce.getLength();

			switch (ce.getOperator()) {
			case D:
				perfectMatch = false;
				DeletionVariation dv = new DeletionVariation();
				dv.setLength(ceLen);
				if (!cramRecord.isNegativeStrand())
					dv.setPosition((int) (readLength - (posInRead)));
				else
					dv.setPosition(posInRead + 1);
				posInSeq += ceLen;
				if (cramRecord.getDeletionVariations() == null)
					cramRecord
							.setDeletionVariations(new LinkedList<DeletionVariation>());
				cramRecord.getDeletionVariations().add(dv);
				break;
			case S:
				posInRead += ceLen;
//				for (int i = 0; i < ceLen; i++) {
//					if (bases[posInRead] != refBases[(int) (cramRecord
//							.getAlignmentStart() + posInSeq) - 1]) {
//						perfectMatch = false;
//						SubstitutionVariation sv = new SubstitutionVariation();
//						if (!cramRecord.isNegativeStrand())
//							sv.setPosition((int) (readLength - (posInRead)));
//						else
//							sv.setPosition(posInRead + 1);
//						sv.setBase(bases[posInRead]);
//						sv.setRefernceBase(refBases[(int) (cramRecord
//								.getAlignmentStart() + posInSeq) - 1]);
//						sv.setQualityScore(qualityScore[posInRead]);
//						if (cramRecord.getSubstitutionVariations() == null)
//							cramRecord
//									.setSubstitutionVariations(new LinkedList<SubstitutionVariation>());
//						cramRecord.getSubstitutionVariations().add(sv);
//					}
//					posInRead++;
//				}
				break;
			case I:
				perfectMatch = false;
				InsertionVariation iv = new InsertionVariation();
				if (!cramRecord.isNegativeStrand())
					iv.setPosition((int) (readLength - (posInRead)));
				else
					iv.setPosition(posInRead + 1);
				iv.setSequence(Arrays.copyOfRange(bases, posInRead, posInRead
						+ ceLen));
				posInRead += ceLen;
				if (cramRecord.getInsertionVariations() == null)
					cramRecord
							.setInsertionVariations(new LinkedList<InsertionVariation>());
				cramRecord.getInsertionVariations().add(iv);
				break;
			case M:
			case X:
				// case S:
				for (int i = 0; i < ceLen; i++) {
					if (bases[posInRead] != refBases[(int) (cramRecord
							.getAlignmentStart() + posInSeq) - 1]) {
						perfectMatch = false;
						SubstitutionVariation sv = new SubstitutionVariation();
						if (!cramRecord.isNegativeStrand())
							sv.setPosition((int) (readLength - (posInRead)));
						else
							sv.setPosition(posInRead + 1);
						sv.setBase(bases[posInRead]);
						sv.setRefernceBase(refBases[(int) (cramRecord
								.getAlignmentStart() + posInSeq) - 1]);
						sv.setQualityScore(qualityScore[posInRead]);
						if (cramRecord.getSubstitutionVariations() == null)
							cramRecord
									.setSubstitutionVariations(new LinkedList<SubstitutionVariation>());
						cramRecord.getSubstitutionVariations().add(sv);
					}
					posInRead++;
					posInSeq++;
				}
				break;
			default:
				break;
			}
		}
		cramRecord.setPerfectMatch(perfectMatch);
	}

	public static CramRecord newRecord2(int alilgnmetStart,
			boolean negativeStrand, List<CigarElement> cigarElements,
			byte[] bases, byte[] qualityScore, byte[] refBases) {
		CramRecord cramRecord = new CramRecord();
		cramRecord.setAlignmentStart(alilgnmetStart);
		cramRecord.setNegativeStrand(negativeStrand);

		addVariations(cramRecord, cigarElements, bases, qualityScore, refBases);
		return cramRecord;
	}

	public static CramRecord newRecord3(SAMRecord record, byte[] refBases) {
		CramRecord cramRecord = new CramRecord();
		cramRecord.setAlignmentStart(record.getAlignmentStart());
		cramRecord.setNegativeStrand(record.getReadNegativeStrandFlag());
		cramRecord.setReadMapped(!record.getReadUnmappedFlag());
		cramRecord.setReadLength(record.getReadLength());
		if (!record.getReadPairedFlag())
			cramRecord.setLastFragment(false);
		else {
			if (record.getFirstOfPairFlag())
				cramRecord.setLastFragment(false);
			else if (record.getSecondOfPairFlag())
				cramRecord.setLastFragment(true);
			else
				throw new IllegalArgumentException(
						"Could not determine if the read is last fragment in the template: "
								+ record.toString());
		}

		addVariations(cramRecord, record.getCigar().getCigarElements(),
				record.getReadBases(), record.getBaseQualities(), refBases);
		return cramRecord;
	}

	@Deprecated
	public static CramRecord newRecord(ReferenceSequenceFile refSeqFile,
			SAMRecord record) {
		CramRecord cramRecord = new CramRecord();
		cramRecord.setAlignmentStart(record.getAlignmentStart());
		cramRecord.setNegativeStrand(record.getReadNegativeStrandFlag());

		boolean perfectMatch = true;
		Cigar cigar = record.getCigar();

		int seqLen = record.getReadLength();

		for (CigarElement ce : cigar.getCigarElements()) {
			switch (ce.getOperator()) {
			case I:
			case S:
				seqLen -= ce.getLength();
				break;
			case D:
				seqLen += ce.getLength();
			default:
				break;
			}
		}
		ReferenceSequence refSeq = refSeqFile.getSubsequenceAt(
				record.getReferenceName(), record.getAlignmentStart(),
				record.getAlignmentStart() + seqLen - 1);
		byte[] refBases = refSeq.getBases();

		int posInRead = 0;
		int posInSeq = 0;
		for (CigarElement ce : cigar.getCigarElements()) {
			switch (ce.getOperator()) {
			case D:
				perfectMatch = false;
				DeletionVariation dv = new DeletionVariation();
				dv.setLength(ce.getLength());
				dv.setPosition(posInRead + 1);
				posInSeq += ce.getLength();
				if (cramRecord.getDeletionVariations() == null)
					cramRecord
							.setDeletionVariations(new LinkedList<DeletionVariation>());
				cramRecord.getDeletionVariations().add(dv);
				break;
			case I:
			case S:
				perfectMatch = false;
				InsertionVariation iv = new InsertionVariation();
				iv.setPosition(posInRead + 1);
				iv.setSequence(Arrays.copyOfRange(record.getReadBases(),
						posInRead, posInRead + ce.getLength()));
				posInRead += ce.getLength();
				if (cramRecord.getInsertionVariations() == null)
					cramRecord
							.setInsertionVariations(new LinkedList<InsertionVariation>());
				cramRecord.getInsertionVariations().add(iv);
				break;
			case M:
			case X:

				for (int i = 0; i < ce.getLength(); i++) {
					if (record.getReadBases()[posInRead] != refBases[posInSeq]) {
						perfectMatch = false;
						SubstitutionVariation sv = new SubstitutionVariation();
						sv.setPosition(posInRead + 1);
						sv.setBase(record.getReadBases()[posInRead]);
						sv.setRefernceBase(refBases[posInSeq]);
						if (cramRecord.getSubstitutionVariations() == null)
							cramRecord
									.setSubstitutionVariations(new LinkedList<SubstitutionVariation>());
						cramRecord.getSubstitutionVariations().add(sv);
					}
					posInRead++;
					posInSeq++;
				}
				break;
			default:
				break;
			}
		}

		cramRecord.setPerfectMatch(perfectMatch);
		return cramRecord;
	}

}
