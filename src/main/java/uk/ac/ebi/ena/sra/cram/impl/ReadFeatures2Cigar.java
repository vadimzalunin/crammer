package uk.ac.ebi.ena.sra.cram.impl;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.InsertBase;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.ReadBase;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;

public class ReadFeatures2Cigar {

	@Deprecated
	public Cigar getCigar(Collection<ReadFeature> features, int readLength) {
		if (features == null || features.isEmpty()) {
			CigarElement ce = new CigarElement(readLength, CigarOperator.M);
			return new Cigar(Arrays.asList(ce));
		}

		List<CigarElement> list = new ArrayList<CigarElement>();
		int prevPos = 1;
		CigarElement ce;
		for (ReadFeature f : features) {
			switch (f.getOperator()) {
			case InsertionVariation.operator:
				InsertionVariation iv = (InsertionVariation) f;

				if (iv.getPosition() - prevPos > 0) {
					ce = new CigarElement(iv.getPosition() - prevPos,
							CigarOperator.M);
					list.add(ce);
				}

				ce = new CigarElement(iv.getSequence().length, CigarOperator.I);
				list.add(ce);
				prevPos = iv.getPosition() + iv.getSequence().length;
				break;
			case InsertBase.operator:
				InsertBase ib = (InsertBase) f;

				if (ib.getPosition() - prevPos > 0) {
					ce = new CigarElement(ib.getPosition() - prevPos,
							CigarOperator.M);
					list.add(ce);
				}

				ce = new CigarElement(1, CigarOperator.I);
				list.add(ce);
				prevPos = ib.getPosition() + 1;
				break;
			case DeletionVariation.operator:
				DeletionVariation dv = (DeletionVariation) f;

				if (dv.getPosition() - prevPos > 0) {
					ce = new CigarElement(dv.getPosition() - prevPos,
							CigarOperator.M);
					list.add(ce);
				}

				ce = new CigarElement(dv.getLength(), CigarOperator.D);
				list.add(ce);
				prevPos = dv.getPosition();
				break;
			default:
				break;
			}
		}

		if (list.isEmpty()) {
			ce = new CigarElement(readLength, CigarOperator.M);
			return new Cigar(Arrays.asList(ce));
		}

		if (readLength - prevPos >= 0) {
			ce = new CigarElement(readLength - prevPos + 1, CigarOperator.M);
			list.add(ce);
		}

		rewriteFlankingInsertsAsSoftClips (list) ;
		return new Cigar(list);
	}
	
	private static final void rewriteFlankingInsertsAsSoftClips(
			List<CigarElement> elements) {
		if (elements.isEmpty())
			return;

		CigarElement first = elements.get(0);
		if (first.getOperator() == CigarOperator.INSERTION) {
			elements.set(0, new CigarElement(first.getLength(),
					CigarOperator.SOFT_CLIP));
		}

		CigarElement last = elements.get(elements.size() - 1);
		if (last.getOperator() == CigarOperator.INSERTION) {
			elements.set(elements.size() - 1, new CigarElement(last.getLength(),
					CigarOperator.SOFT_CLIP));
		}
	}

	private static final CigarOperator detectCigarOperator(ReadFeature f) {
		switch (f.getOperator()) {
		case InsertionVariation.operator:
		case InsertBase.operator:
			return CigarOperator.INSERTION;
		case DeletionVariation.operator:
			return CigarOperator.DELETION;
		default:
			return CigarOperator.MATCH_OR_MISMATCH;
		}
	}

	public Cigar getCigar2(Collection<ReadFeature> features, int readLength) {
		if (features == null || features.isEmpty()) {
			CigarElement ce = new CigarElement(readLength, CigarOperator.M);
			return new Cigar(Arrays.asList(ce));
		}

		List<CigarElement> list = new ArrayList<CigarElement>();
		int totalOpLen = 1;
		CigarElement ce;
		CigarOperator lastOperator = CigarOperator.MATCH_OR_MISMATCH;
		int lastOpLen = 0;
		int lastOpPos = 1;
		CigarOperator co = null;
		int rfLen = 0;
		for (ReadFeature f : features) {
			// if (lastOperator == CigarOperator.DELETION)

			int gap = f.getPosition() - (lastOpPos + lastOpLen);
			if (gap > 0) {
				if (lastOperator != CigarOperator.MATCH_OR_MISMATCH) {
					list.add(new CigarElement(lastOpLen, lastOperator));
					lastOpPos += lastOpLen;
					totalOpLen += lastOpLen;
					lastOpLen = gap;
				} else {
					lastOpLen += gap;
				}

				lastOperator = CigarOperator.MATCH_OR_MISMATCH;
			}

			switch (f.getOperator()) {
			case InsertionVariation.operator:
				co = CigarOperator.INSERTION;
				rfLen = ((InsertionVariation) f).getSequence().length;
				break;
			case InsertBase.operator:
				co = CigarOperator.INSERTION;
				rfLen = 1;
				break;
			case DeletionVariation.operator:
				co = CigarOperator.DELETION;
				rfLen = ((DeletionVariation) f).getLength();
				break;
			case SubstitutionVariation.operator:
			case ReadBase.operator:
				co = CigarOperator.MATCH_OR_MISMATCH;
				rfLen = 1;
				break;
			default:
				continue;
			}

			if (lastOperator != co) {
				// add last feature
				if (lastOpLen > 0) {
					list.add(new CigarElement(lastOpLen, lastOperator));
					totalOpLen += lastOpLen;
				}
				lastOperator = co;
				lastOpLen = rfLen;
				lastOpPos = f.getPosition();
				if (co == CigarOperator.DELETION)
					lastOpPos -= rfLen;
			} else
				lastOpLen += rfLen;
		}

		if (lastOperator != null) {
			if (lastOperator != CigarOperator.M) {
				list.add(new CigarElement(lastOpLen, lastOperator));
				if (readLength >= lastOpPos + lastOpLen) {
					ce = new CigarElement(readLength - (lastOpLen + lastOpPos)+1,
							CigarOperator.M);
					list.add(ce);
				}
			} else if (readLength > lastOpPos - 1) {
				ce = new CigarElement(readLength - lastOpPos + 1,
						CigarOperator.M);
				list.add(ce);
			}
		}

		if (list.isEmpty()) {
			ce = new CigarElement(readLength, CigarOperator.M);
			return new Cigar(Arrays.asList(ce));
		}

		rewriteFlankingInsertsAsSoftClips (list) ;
		return new Cigar(list);
	}

}
