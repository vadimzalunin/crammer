package uk.ac.ebi.ena.sra.cram.impl;

import java.util.Collections;

import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;

public class ReverseVariations {

	public void reverse(CramRecord record) {
		if (record.getReadFeatures() == null
				|| record.getReadFeatures().isEmpty())
			return;

		for (ReadFeature f : record.getReadFeatures())
			f.setPosition((int) (record.getReadLength() - f.getPosition() + 1));

		Collections.reverse(record.getReadFeatures());
	}
}
