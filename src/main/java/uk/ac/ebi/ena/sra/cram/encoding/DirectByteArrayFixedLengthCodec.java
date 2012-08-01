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

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class DirectByteArrayFixedLengthCodec implements BitCodec<byte[]> {
	private int length;

	public DirectByteArrayFixedLengthCodec(int length) {
		this.length = length;
	}

	@Override
	public byte[] read(BitInputStream bis) throws IOException {
		byte[] bytes = new byte[length];
		bis.readAlignedBytes(bytes);

		return bytes;
	}

	@Override
	public long write(BitOutputStream bos, byte[] object) throws IOException {
		if (object.length != length)
			throw new RuntimeException("Array length must be " + length);
		bos.write(object);
		return 8 * length;
	}

	@Override
	public long numberOfBits(byte[] object) {
		return 8 * length;
	}

}
