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

public class ByteArrayHuffmanCodec implements BitCodec<byte[]> {
	private HuffmanByteCodec2 byteCodec;
	private ByteBuffer buf = ByteBuffer.allocate(1000);
	private byte stopByte;

	public ByteArrayHuffmanCodec(byte[] alphabet, int[] freqs, byte stopByte) {
		this.stopByte = stopByte;
		HuffmanTree<Byte> tree = HuffmanCode.buildTree(freqs, Utils.autobox(alphabet));
		byteCodec = new HuffmanByteCodec2(tree);
	}

	@Override
	public byte[] read(BitInputStream bis) throws IOException {
		buf.clear();
		byte b;
		while ((b = byteCodec.read(bis)) != stopByte)
			buf.put(b);

		byte[] sequence = new byte[buf.position()];
		buf.flip();
		buf.get(sequence);
		return sequence;
	}

	@Override
	public long write(BitOutputStream bos, byte[] bytes) throws IOException {
		if (bytes.length == 0)
			return 0;

		long len = 0;
		for (byte b : bytes)
			len += byteCodec.write(bos, b);

		if (bytes[bytes.length - 1] != stopByte)
			len += byteCodec.write(bos, stopByte);

		return len;
	}

	@Override
	public long numberOfBits(byte[] bytes) {
		if (bytes.length == 0)
			return 0;

		long len = 0;
		for (byte b : bytes)
			len += byteCodec.numberOfBits(b);

		if (bytes[bytes.length - 1] != stopByte)
			len += byteCodec.numberOfBits(stopByte);

		return len;
	}

}
