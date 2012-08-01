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

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class VariableLengthByteArrayHuffmanCodec implements BitCodec<byte[]> {
	private HuffmanByteCodec2 byteCodec;
	private HuffmanCodec<Integer> lenCodec;

	public VariableLengthByteArrayHuffmanCodec(byte[] alphabet, int[] freqs, int[] lenAlphabet, int[] lenFreqs) {

		int maxLen = 0;
		for (int len : lenAlphabet)
			if (maxLen < len)
				maxLen = len;

		HuffmanTree<Byte> tree = HuffmanCode.buildTree(freqs, Utils.autobox(alphabet));
		byteCodec = new HuffmanByteCodec2(tree);

		HuffmanTree<Integer> lenTree = HuffmanCode.buildTree(lenFreqs, Utils.autobox(lenAlphabet));
		lenCodec = new HuffmanCodec<Integer>(lenTree);
	}

	@Override
	public byte[] read(BitInputStream bis) throws IOException {
		Integer len = lenCodec.read(bis) ;
		byte[] sequence = new byte[len];

		for (int i = 0; i < sequence.length; i++)
			sequence[i] = byteCodec.read(bis);

		return sequence;
	}

	@Override
	public long write(BitOutputStream bos, byte[] bytes) throws IOException {
		long len = 0;

		len += lenCodec.write(bos, bytes.length);

		for (byte b : bytes)
			len += byteCodec.write(bos, b);

		return len;
	}

	@Override
	public long numberOfBits(byte[] bytes) {
		long len = 0;

		len += lenCodec.numberOfBits(bytes.length);

		for (byte b : bytes)
			len += byteCodec.numberOfBits(b);

		return len;
	}

}
