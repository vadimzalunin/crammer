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

import uk.ac.ebi.ena.sra.cram.encoding.GammaCodec;

class GammaCodecStub extends GammaCodec implements NumberCodecStub {

	public GammaCodecStub() {
		super();
	}

	@Override
	public EncodingAlgorithm getEncoding() {
		return EncodingAlgorithm.GAMMA;
	}

	@Override
	public String getStringRepresentation() {
		return String.format("%d,%d", getOffset(), isLenCodingBit() ? 1 : 0);
	}

	@Override
	public void initFromString(String spec) throws CramCompressionException {
		String[] params = StringRepresentation.parse(spec);
		switch (params.length) {
		case 0:
			break;
		case 2:
			long offset = StringRepresentation.toLong(params[0]);
			setOffset(offset);
			boolean lenCodingBit = StringRepresentation.toBoolean(params[1]);
			setLenCodingBit(lenCodingBit);
			break;
		default:
			throw new CramCompressionException(
					"Not supported number of parameters to unary codec: "
							+ params.length);
		}
	}

	@Override
	public Object[] getParameters() {
		return new Object[] { getOffset(), isLenCodingBit() };
	}

	@Override
	public void setParameters(Object[] params) {
		setOffset((Long) params[0]);
		setLenCodingBit((Boolean) params[1]);
	}
}
