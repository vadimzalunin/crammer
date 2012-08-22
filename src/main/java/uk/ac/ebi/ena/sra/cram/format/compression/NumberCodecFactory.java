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

public class NumberCodecFactory {

	public static NumberCodecStub createStub(EncodingAlgorithm encoding) throws CramCompressionException {
		switch (encoding) {
		case NULL:
			return new NullCodecStub();
		case GOLOMB:
			return new GolombCodecStub();
		case UNARY:
			return new UnaryCodecStub();
		case GAMMA:
			return new GammaCodecStub();
		case GOLOMB_RICE:
			return new GolombRiceCodecStub();
		case SUBEXP:
			return new SubexpCodecStub();
		case BETA:
			return new BetaCodecStub();

		default:
			break;
		}
		throw new CramCompressionException("Unsupported codec: " + encoding);
	}
}
