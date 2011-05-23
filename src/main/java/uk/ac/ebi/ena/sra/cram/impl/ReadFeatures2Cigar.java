package uk.ac.ebi.ena.sra.cram.impl;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;

public class ReadFeatures2Cigar {

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

		return new Cigar(list);
	}
}
