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

public class ByteArraySequenceBaseProvider implements SequenceBaseProvider {
	private byte[] sequence;
	private boolean circular = false;
	private int N_extension = 0;

	public ByteArraySequenceBaseProvider(byte[] sequence) {
		this.sequence = sequence;
	}

	@Override
	public byte getBaseAt(String seqName, long position) {
		if (position < sequence.length)
			return sequence[(int) position];

		if (circular)
			return sequence[(int) (position % sequence.length)];

		if (position < sequence.length + N_extension)
			return 'N';

		throw new RuntimeException(String.format(
				"Reference position out of range: in sequence %s, length %d, position %d.", seqName, sequence.length,
				position));
	}

	@Override
	public void copyBases(String sequenceName, long from, int len, byte[] dest) throws IOException {
		try {
			if (from + len > sequence.length) {
				Arrays.fill(dest, (byte) 'N');
				System.arraycopy(sequence, (int) from, dest, 0, sequence.length - (int) from);
			} else
				System.arraycopy(sequence, (int) from, dest, 0, len);
		} catch (ArrayIndexOutOfBoundsException e) {
			System.err.printf("Offensive request: sequence %s from=%d len=%d\n", sequenceName, from, len);
			throw e;
		}
	}

	public boolean isCircular() {
		return circular;
	}

	public void setCircular(boolean circular) {
		this.circular = circular;
	}

	public int getN_extension() {
		return N_extension;
	}

	public void setN_extension(int n_extension) {
		N_extension = n_extension;
	}

}
