package uk.ac.ebi.ena.sra.cram.impl;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;

public class RestoreQualityScores {

	private byte defaultQualityScore = 33 + 30;

	public byte[] restoreQualityScores(CramRecord record) throws IOException {
		byte[] scores = new byte[(int) record.getReadLength()];
		Arrays.fill(scores, defaultQualityScore);

		if (record.getReadFeatures() == null
				|| record.getReadFeatures().isEmpty())
			return scores;

		List<ReadFeature> variations = record.getReadFeatures();
		for (ReadFeature v : variations) {
			if (v.getOperator() == SubstitutionVariation.operator) {
				SubstitutionVariation sv = (SubstitutionVariation) v;
				scores[sv.getPosition() - 1] = sv.getQualityScore();
			}
		}
		return scores;
	}
}
