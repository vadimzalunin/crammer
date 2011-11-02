package uk.ac.ebi.ena.sra.cram.impl;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import uk.ac.ebi.ena.sra.cram.format.BaseQualityScore;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.ReadBase;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;

public class RestoreQualityScores {

	private byte defaultQualityScore = 32;

	public byte[] restoreQualityScores(CramRecord record) throws IOException {
		byte[] scores = new byte[(int) record.getReadLength()];
		Arrays.fill(scores, defaultQualityScore);

		if (record.getReadFeatures() == null
				|| record.getReadFeatures().isEmpty())
			return scores;

		List<ReadFeature> variations = record.getReadFeatures();
		for (ReadFeature v : variations) {
			switch (v.getOperator()) {
			case BaseQualityScore.operator:
				BaseQualityScore bqs = (BaseQualityScore) v ;
				scores[bqs.getPosition() - 1] = bqs.getQualityScore() ;
				break ;
//			case SubstitutionVariation.operator:
//				SubstitutionVariation sv = (SubstitutionVariation) v;
//				scores[sv.getPosition() - 1] = sv.getQualityScore();
//				break;

			default:
				break;
			}
		}

		// ReadBase has more weight: 
		for (ReadFeature v : variations) {
			switch (v.getOperator()) {
			case ReadBase.operator:
				ReadBase base = (ReadBase) v;
				scores[base.getPosition() - 1] = base.getQualityScore();
				break;

			default:
				break;
			}
		}
		
		record.setQualityScores(scores) ;
		return scores;
	}
}
