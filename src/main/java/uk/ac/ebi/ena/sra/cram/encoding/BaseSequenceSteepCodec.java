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

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

class BaseSequenceSteepCodec implements BitCodec<byte[]> {
	private byte[] order = "ACGTNS".getBytes();
	private BitCode[] codes;
	private ByteBuffer buffer = ByteBuffer.allocate(1024);

	public BaseSequenceSteepCodec(byte[] order) {
		if (order.length != 6)
			throw new IllegalArgumentException("Expecting 5 bases and 1 stop only but got: " + new String(order));

		this.order = order;
		this.codes = new BitCode[255];

		codes[order[0]] = new BitCode(0, 1);
		codes[order[1]] = new BitCode(2, 2);
		codes[order[2]] = new BitCode(6, 3);
		codes[order[3]] = new BitCode(14, 4);
		codes[order[4]] = new BitCode(30, 5);
		codes[order[5]] = new BitCode(31, 5);
	}

	@Override
	public byte[] read(BitInputStream bis) throws IOException {
		buffer.clear();

		while (true) {
			int trueBitCounter = 0;
			int bitsReadCounter = 0;

			while (bitsReadCounter++ < order.length - 1 && bis.readBit())
				trueBitCounter++;

			if (order[trueBitCounter] == 'S')
				break;
			buffer.put(order[trueBitCounter]);
		}

		buffer.flip();
		byte[] seq = new byte[buffer.limit()];
		buffer.get(seq);
		return seq;
	}

	@Override
	public long write(BitOutputStream bos, byte[] bases) throws IOException {
		int length = 0;
		for (byte base : bases) {
			if (base == 'S')
				throw new IllegalArgumentException("Stop code is not allowed in a base sequence.");
			BitCode code = codes[base];
			bos.write(code.value, code.bits);
			length += code.bits;
		}

		bos.write(codes['S'].value, codes['S'].bits);
		length += codes['S'].bits;

		return length;
	}

	@Override
	public long numberOfBits(byte[] bases) {
		int length = 0;
		for (byte base : bases) {
			if (base == 'S')
				throw new IllegalArgumentException("Stop code is not allowed in a base sequence.");
			BitCode code = codes[base];
			length += code.bits;
		}

		length += codes['S'].bits;

		return length;
	}
}
