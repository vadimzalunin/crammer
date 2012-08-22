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

public class HuffmanCodec<V> implements BitCodec<V> {
	private HuffmanTree<V> tree;
	private TreeMap<V, HuffmanBitCode<V>> codes;

	public HuffmanCodec(HuffmanTree<V> tree) {
		super();
		this.tree = tree;
		codes = new TreeMap<V, HuffmanBitCode<V>>();
		getBitCode(tree, new HuffmanBitCode<V>(), codes);
		// System.out.println("Dumping Huffman codec: ");
		// for (V v : codes.keySet())
		// System.out.println(v + ": " + codes.get(v));
		// for (byte i = 0; i < 100; i++)
		// System.out.println(i + ": " + codes.get(i));
	}

	@Override
	public V read(BitInputStream bis) throws IOException {
		HuffmanLeaf<V> leaf = (HuffmanLeaf<V>) findLeaf(bis, tree);
		return leaf.value;
	}

	private HuffmanTree<V> findLeaf(BitInputStream bis, HuffmanTree<V> tree) throws IOException {
		if (tree instanceof HuffmanLeaf)
			return tree;

		HuffmanNode<V> node = (HuffmanNode<V>) tree;
		return findLeaf(bis, bis.readBit() ? node.right : node.left);
	}

	@Override
	public long write(BitOutputStream bos, V symbol) throws IOException {
		HuffmanBitCode<V> bitCode = codes.get(symbol);
		if (bitCode == null)
			throw new RuntimeException("Huffman code not found for value: " + symbol);
		bos.write(bitCode.bitCode, bitCode.bitLentgh);
		return bitCode.bitLentgh;
	}

	@Override
	public long numberOfBits(V symbol) {
		return codes.get(symbol).bitLentgh;
	}

	private static class HuffmanBitCode<T> {
		long bitCode;
		int bitLentgh;
	}

	private static <T> void getBitCode(HuffmanTree<T> tree, HuffmanBitCode<T> code, Map<T, HuffmanBitCode<T>> codes) {
		if (tree instanceof HuffmanLeaf) {
			HuffmanLeaf<T> leaf = (HuffmanLeaf<T>) tree;
			HuffmanBitCode<T> readyCode = new HuffmanBitCode<T>();
			readyCode.bitCode = code.bitCode;
			readyCode.bitLentgh = code.bitLentgh;
			codes.put(leaf.value, readyCode);
			return;

		} else if (tree instanceof HuffmanNode) {
			HuffmanNode<T> node = (HuffmanNode<T>) tree;

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

	public HuffmanTree<V> getTree() {
		return tree;
	}

	public void clear() {
		codes.clear();
		codes = null;
		tree = null;
	}
}
