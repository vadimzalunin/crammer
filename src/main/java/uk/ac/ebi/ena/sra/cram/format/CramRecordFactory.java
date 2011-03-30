package uk.ac.ebi.ena.sra.cram.format;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;

public class CramRecordFactory {

	private static void addVariations(CramRecord cramRecord,
			List<CigarElement> cigarElements, byte[] bases, byte[] refBases) {
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
				dv.setPosition(posInRead + 1);
				posInSeq += ceLen;
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

				for (int i = 0; i < ceLen; i++) {
					if (bases[posInRead] != refBases[(int) (cramRecord
							.getAlignmentStart() + posInSeq)]) {
						perfectMatch = false;
						SubstitutionVariation sv = new SubstitutionVariation();
						sv.setPosition(posInRead + 1);
						sv.setBase(bases[posInRead]);
						sv.setRefernceBase(refBases[(int) (cramRecord
								.getAlignmentStart() + posInSeq)]);
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
			byte[] bases, byte[] refBases) {
		CramRecord cramRecord = new CramRecord();
		cramRecord.setAlignmentStart(alilgnmetStart);
		cramRecord.setNegativeStrand(negativeStrand);

		addVariations(cramRecord, cigarElements, bases, refBases);
		return cramRecord;
	}

	public static CramRecord newRecord3(SAMRecord record, byte[] refBases) {
		CramRecord cramRecord = new CramRecord();
		cramRecord.setAlignmentStart(record.getAlignmentStart());
		cramRecord.setNegativeStrand(record.getReadNegativeStrandFlag());

		addVariations(cramRecord, record.getCigar().getCigarElements(),
				record.getReadBases(), refBases);
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
