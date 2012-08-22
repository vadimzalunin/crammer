/*******************************************************************************
 * Copyright 2012 EMBL-EBI, Hinxton outstation
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *   http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/
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
		if (record.getQualityScores() != null)
			return record.getQualityScores();

		byte[] scores = new byte[(int) record.getReadLength()];
		Arrays.fill(scores, defaultQualityScore);

		if (record.getReadFeatures() == null || record.getReadFeatures().isEmpty())
			return scores;

		List<ReadFeature> variations = record.getReadFeatures();
		for (ReadFeature v : variations) {
			switch (v.getOperator()) {
			case BaseQualityScore.operator:
				BaseQualityScore bqs = (BaseQualityScore) v;
				scores[bqs.getPosition() - 1] = bqs.getQualityScore();
				break;
			// case SubstitutionVariation.operator:
			// SubstitutionVariation sv = (SubstitutionVariation) v;
			// scores[sv.getPosition() - 1] = sv.getQualityScore();
			// break;

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

		record.setQualityScores(scores);
		return scores;
	}
}
