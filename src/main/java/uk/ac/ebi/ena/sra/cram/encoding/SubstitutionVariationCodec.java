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
package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.NullBitOutputStream;

public class SubstitutionVariationCodec implements BitCodec<SubstitutionVariation> {
	public BitCodec<BaseChange> baseChangeCodec;

	// public BitCodec<Byte> qualityScoreCodec;

	@Override
	public SubstitutionVariation read(BitInputStream bis) throws IOException {
		SubstitutionVariation v = new SubstitutionVariation();
		// position is not read here because we need to keep track of previous
		// values read from the codec. See ReadFeatureCodec.
		long position = -1L;
		BaseChange baseChange = baseChangeCodec.read(bis);
		// byte qualityScore = qualityScoreCodec.read(bis);

		v.setPosition((int) position);
		v.setBaseChange(baseChange);
		// v.setQualityScore(qualityScore);

		return v;
	}

	@Override
	public long write(BitOutputStream bos, SubstitutionVariation v) throws IOException {
		long len = 0L;

		BaseChange baseChange = v.getBaseChange();
		if (baseChange == null)
			baseChange = new BaseChange(v.getRefernceBase(), v.getBase());

		long baseChangeLen = 0L;
		baseChangeLen -= len;
		len += baseChangeCodec.write(bos, baseChange);
		baseChangeLen += len;

		// len += qualityScoreCodec.write(bos, v.getQualityScore());

		return len;
	}

	@Override
	public long numberOfBits(SubstitutionVariation v) {
		try {
			return write(NullBitOutputStream.INSTANCE, v);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

}
