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
import java.nio.ByteBuffer;
import java.util.Arrays;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

class BaseSequenceFlatCodec implements BitCodec<byte[]> {
	private byte[] order = "ACGTNS".getBytes();
	private int[] base2indexArray;
	private ByteBuffer buffer = ByteBuffer.allocate(1024);

	public BaseSequenceFlatCodec(byte[] order) {
		if (order.length != 6)
			throw new IllegalArgumentException(
					"Expecting 5 bases and 1 stop only but got: "
							+ new String(order));

		this.order = order;
		this.base2indexArray = new int[255];
		Arrays.fill(base2indexArray, -1);
		for (int i = 0; i < 255; i++)
			base2indexArray[i] = -1;

		int index = 0;
		for (byte base : order)
			base2indexArray[base] = index++;
	}

	@Override
	public byte[] read(BitInputStream bis) throws IOException {
		buffer.clear();

		while (true) {
			int threeBits = bis.readBits(3);

			if (threeBits >= order.length)
				throw new RuntimeException("Unexpected base flat code: "
						+ threeBits);

			if (order[threeBits] == 'S')
				break;
			buffer.put(order[threeBits]);
		}

		buffer.flip() ;
		byte[] seq = new byte[buffer.limit()];
		buffer.get(seq);
		return seq;
	}

	@Override
	public long write(BitOutputStream bis, byte[] bases) throws IOException {
		for (byte base : bases)
			bis.write(base2indexArray[base], 3);

		bis.write(base2indexArray['S'], 3);
		return (bases.length + 1) * 3;
	}

	@Override
	public long numberOfBits(byte[] bases) {
		return (bases.length + 1) * 3;
	}
}
