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

import uk.ac.ebi.ena.sra.cram.encoding.UnaryCodec;

class UnaryCodecStub extends UnaryCodec implements NumberCodecStub {

	public UnaryCodecStub() {
		super();
	}

	@Override
	public EncodingAlgorithm getEncoding() {
		return EncodingAlgorithm.UNARY;
	}

	@Override
	public String getStringRepresentation() {
		return String.format("%d,%d", getOffset(), isStopBit() ? 1 : 0);
	}

	@Override
	public void initFromString(String spec) throws CramCompressionException {
		String[] params = StringRepresentation.parse(spec);
		switch (params.length) {
		case 0:
			break;
		case 2:
			boolean stopBit = StringRepresentation.toBoolean(params[0]);
			setStopBit(stopBit);
			long offset = StringRepresentation.toLong(params[1]);
			setOffset(offset);
			break;
		default:
			throw new CramCompressionException("Not supported number of parameters to unary codec: " + params.length);
		}
	}

	@Override
	public Object[] getParameters() {
		return new Object[] { getOffset(), isStopBit() };
	}

	@Override
	public void setParameters(Object[] params) {
		setOffset((Long) params[0]);
		setStopBit((Boolean) params[1]);
	}
}
