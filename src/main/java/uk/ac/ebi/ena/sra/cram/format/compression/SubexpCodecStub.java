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

import uk.ac.ebi.ena.sra.cram.encoding.SubexpCodec;

class SubexpCodecStub extends SubexpCodec implements NumberCodecStub {

	public SubexpCodecStub() {
		super(1);
	}

	@Override
	public EncodingAlgorithm getEncoding() {
		return EncodingAlgorithm.SUBEXP;
	}

	@Override
	public String getStringRepresentation() {
		return String.format("%d,%d,%d", getK(), getOffset(), isUnaryBit() ? 1 : 0);
	}

	@Override
	public void initFromString(String spec) throws CramCompressionException {
		String[] params = StringRepresentation.parse(spec);
		switch (params.length) {
		case 1:
			int k = StringRepresentation.toInt(params[0]);
			setK(k);
			break;
		case 3:
			k = StringRepresentation.toInt(params[0]);
			setK(k);
			long offset = StringRepresentation.toLong(params[1]);
			setOffset(offset);
			boolean quotientBit = StringRepresentation.toBoolean(params[2]);
			setUnaryBit(quotientBit);
			break;
		default:
			throw new CramCompressionException("Not supported number of parameters to golomb-rice codec: "
					+ params.length);
		}
	}

	@Override
	public Object[] getParameters() {
		return new Object[] { getK(), getOffset(), isUnaryBit() };
	}

	@Override
	public void setParameters(Object[] params) {
		setK((Long) params[0]);
		setOffset((Long) params[1]);
		setUnaryBit((Boolean) params[2]);
	}
}
