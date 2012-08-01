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
package uk.ac.ebi.ena.sra.cram.format.compression;

import uk.ac.ebi.ena.sra.cram.encoding.BetaCodec;

class NullCodecStub extends BetaCodec implements NumberCodecStub {

	@Override
	public EncodingAlgorithm getEncoding() {
		return EncodingAlgorithm.NULL;
	}

	@Override
	public String getStringRepresentation() {
		return "";
	}

	@Override
	public void initFromString(String spec) throws CramCompressionException {
		if (spec != null && spec.length() > 1)
			throw new CramCompressionException(
					"Null codec does not expect any parameteres: " + spec);
	}

	@Override
	public Object[] getParameters() {
		return new Object[] {};
	}

	@Override
	public void setParameters(Object[] params) {
	}
}
