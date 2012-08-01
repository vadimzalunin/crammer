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
import java.util.Map;
import java.util.TreeMap;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanLeaf;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanNode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class HuffmanByteCodec implements BitCodec<Byte> {
	private HuffmanTree<Byte> tree;
	private TreeMap<Byte, HuffmanBitCode> codes;
	private HuffmanBitCode[] bitCodes = new HuffmanBitCode[256] ;

	public HuffmanByteCodec(HuffmanTree<Byte> tree) {
		super();
		this.tree = tree;
		codes = new TreeMap<Byte, HuffmanBitCode>();
		getBitCode(tree, new HuffmanBitCode(), codes);
		
		for (Map.Entry<Byte, HuffmanBitCode> entry: codes.entrySet()) 
			bitCodes[entry.getKey() & 0xFF] = entry.getValue() ;
	}

	@Override
	public Byte read(BitInputStream bis) throws IOException {
		HuffmanLeaf<Byte> leaf = (HuffmanLeaf<Byte>) findLeaf(bis, tree);
		return leaf.value;
	}

	private HuffmanTree<Byte> findLeaf(BitInputStream bis, HuffmanTree<Byte> tree)
			throws IOException {
		if (tree instanceof HuffmanLeaf)
			return tree;

		HuffmanNode<Byte> node = (HuffmanNode<Byte>) tree;
		return findLeaf(bis, bis.readBit() ? node.right : node.left);
	}

	@Override
	public long write(BitOutputStream bos, Byte symbol) throws IOException {
		HuffmanBitCode bitCode = bitCodes[symbol & 0xFF] ;
		if (bitCode == null)
			throw new RuntimeException("Huffman code not found for value: "
					+ symbol);
		bos.write(bitCode.bitCode, bitCode.bitLentgh);
		return bitCode.bitLentgh;
	}

	@Override
	public long numberOfBits(Byte symbol) {
		return codes.get(symbol).bitLentgh;
	}

	private static class HuffmanBitCode {
		int bitCode;
		int bitLentgh;
	}

	private static <T> void getBitCode(HuffmanTree<Byte> tree,
			HuffmanBitCode code, Map<Byte, HuffmanBitCode> codes) {
		if (tree instanceof HuffmanLeaf) {
			HuffmanLeaf<Byte> leaf = (HuffmanLeaf<Byte>) tree;
			HuffmanBitCode readyCode = new HuffmanBitCode();
			readyCode.bitCode = code.bitCode;
			readyCode.bitLentgh = code.bitLentgh;
			codes.put(leaf.value, readyCode);
			return;

		} else if (tree instanceof HuffmanNode) {
			HuffmanNode<Byte> node = (HuffmanNode<Byte>) tree;

			// traverse left
			code.bitCode = code.bitCode << 1;
			code.bitLentgh++;

			getBitCode(node.left, code, codes);
			code.bitCode = code.bitCode >>> 1;
			code.bitLentgh--;

			// traverse right
			code.bitCode = code.bitCode << 1 | 1;
			code.bitLentgh++;

			getBitCode(node.right, code, codes);
			code.bitCode = code.bitCode >>> 1;
			code.bitLentgh--;
		}
	}

	public HuffmanTree<Byte> getTree() {
		return tree;
	}
	
	public void clear () {
		codes.clear() ;
		codes = null ;
		tree = null ;
	}
}
