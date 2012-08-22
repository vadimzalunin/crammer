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

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;

public class ReadFeatureOperatorCodec extends HuffmanByteCodec2 {

	public ReadFeatureOperatorCodec() {
		this(new int[] { 100, 50, 10, 5 }, new Byte[] { 'N', 'S', 'I', 'D' });
	}

	public ReadFeatureOperatorCodec(int[] frequencies, Byte[] alphabet) {
		super(HuffmanCode.buildTree(frequencies, alphabet));
	}

}
