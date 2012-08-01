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

import uk.ac.ebi.ena.sra.cram.encoding.GolombCodec;

class GolombCodecStub extends GolombCodec implements NumberCodecStub {

	public GolombCodecStub() {
		super (2) ;
	}

	@Override
	public EncodingAlgorithm getEncoding() {
		return EncodingAlgorithm.GOLOMB;
	}

	@Override
	public String getStringRepresentation() {
		return String.format("%d,%d,%d", getM(), getOffset(),
				isQuotientBit() ? 1 : 0);
	}

	@Override
	public void initFromString(String spec) throws CramCompressionException {
		String[] params = StringRepresentation.parse(spec);
		switch (params.length) {
		case 1:
			int m = StringRepresentation.toInt(params[0]);
			setM(m);
			break ;
		case 3:
			m = StringRepresentation.toInt(params[0]);
			setM(m);
			long offset = StringRepresentation.toLong(params[1]);
			setOffset(offset);
			boolean quotientBit = StringRepresentation.toBoolean(params[2]);
			setQuotientBit(quotientBit);
			break ;
		default:
			throw new CramCompressionException(
					"Not supported number of parameters to golomb codec: "
							+ params.length);
		}
	}
	
	@Override
	public Object[] getParameters() {
		return new Object[]{getM(), getOffset(), isQuotientBit()};
	}
	@Override
	public void setParameters(Object[] params) {
		setM((Long) params[0]);
		setOffset((Long) params[1]);
		setQuotientBit((Boolean) params[2]);
	}
}
