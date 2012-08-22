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

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class FixedLengthByteArrayHuffmanCodec implements BitCodec<byte[]> {
	private HuffmanByteCodec2 byteCodec;
	private ByteBuffer buf;
	private int length;

	public FixedLengthByteArrayHuffmanCodec(byte[] alphabet, int[] freqs, int len) {
		this.length = len;
		buf = ByteBuffer.allocate(len);
		HuffmanTree<Byte> tree = HuffmanCode.buildTree(freqs, Utils.autobox(alphabet));
		byteCodec = new HuffmanByteCodec2(tree);
	}

	@Override
	public byte[] read(BitInputStream bis) throws IOException {
		buf.clear();
		for (int i = 0; i < length; i++)
			buf.put(byteCodec.read(bis));

		byte[] sequence = new byte[buf.position()];
		buf.flip();
		buf.get(sequence);
		return sequence;
	}

	@Override
	public long write(BitOutputStream bos, byte[] bytes) throws IOException {
		if (bytes.length != length)
			throw new RuntimeException("Number of bytes in the value is different from " + length);

		long len = 0;
		for (byte b : bytes)
			len += byteCodec.write(bos, b);

		return len;
	}

	@Override
	public long numberOfBits(byte[] bytes) {
		if (bytes.length != length)
			throw new RuntimeException("Number of bytes in the value is different from " + length);

		long bitCount = 0;
		for (byte b : bytes)
			bitCount += byteCodec.numberOfBits(b);

		return bitCount;
	}

}
