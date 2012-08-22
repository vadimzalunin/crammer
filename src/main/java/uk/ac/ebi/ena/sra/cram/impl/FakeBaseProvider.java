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

import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;

public class FakeBaseProvider implements SequenceBaseProvider {
	private byte base = 'N';

	@Override
	public byte getBaseAt(String seqName, long position) {
		return base;
	}

	@Override
	public void copyBases(String sequenceName, long from, int len, byte[] dest) throws IOException {
		Arrays.fill(dest, 0, len, base);
	}

}
