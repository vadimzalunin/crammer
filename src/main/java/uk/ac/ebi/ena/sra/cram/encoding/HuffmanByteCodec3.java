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

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanLeaf;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanNode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.format.ByteFrequencies;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

public class HuffmanByteCodec3 implements BitCodec<Byte> {
	private HuffmanTree<Byte> tree;
	private TreeMap<Byte, HuffmanBitCode> codes;
	private HuffmanBitCode[] bitCodes = new HuffmanBitCode[256];
	// assuming code length cannot be more than 1024:
	private Integer[] codeLentghSorted;
	private Map<Integer, Map<Long, Byte>> codeCache = new HashMap<Integer, Map<Long, Byte>>();
	private Map<Long, Byte>[] codeMaps;

	Node root;
	Node[][] nodes;
	private List<Node> nodeList;

	class Node {
		int len = 0;
		HuffmanBitCode[] leaf;
		Node[] children;
	}

	public static HuffmanByteCodec3 build(ByteFrequencies bf) {
		HuffmanTree<Byte> tree = HuffmanCode.buildTree(bf.getFrequencies(), Utils.autobox(bf.getValues()));
		return new HuffmanByteCodec3(tree);
	}

	public HuffmanByteCodec3(HuffmanTree<Byte> tree) {
		super();
		this.tree = tree;
		codes = new TreeMap<Byte, HuffmanBitCode>();
		getBitCode(tree, new HuffmanBitCode(), codes);

		if (codes.isEmpty())
			return;

		for (Map.Entry<Byte, HuffmanBitCode> entry : codes.entrySet())
			bitCodes[entry.getKey() & 0xFF] = entry.getValue();

		Set<Integer> lenMap = new HashSet<Integer>();
		for (int i = 0; i < bitCodes.length; i++) {
			HuffmanBitCode code = bitCodes[i];
			if (code == null)
				continue;
			lenMap.add(code.bitLentgh);
			Map<Long, Byte> codeMap = codeCache.get(code.bitLentgh);
			if (codeMap == null) {
				codeMap = new HashMap<Long, Byte>();
				codeCache.put(code.bitLentgh, codeMap);
			}
			codeMap.put(new Long(code.bitCode), new Byte((byte) i));
		}

		codeLentghSorted = new Integer[lenMap.size()];
		lenMap.toArray(codeLentghSorted);
		Arrays.sort(codeLentghSorted);

		codeMaps = new Map[codeLentghSorted[codeLentghSorted.length - 1] + 1];
		for (int len : codeLentghSorted) {
			codeMaps[len] = codeCache.get(len);
		}

		root = bn(0, codes.values());

		// nodeList = new ArrayList<Node>(lenMap.size());
		// int prevLen = 0;
		// for (int len : codeLentghSorted) {
		// System.out.printf("For length=%d\n", len);
		//
		// Map<Long, Byte> codeMap = codeMaps[len];
		// Node node = new Node();
		// nodeList.add(node);
		//
		// List<HuffmanBitCode> codesOfTheLength = new
		// ArrayList<HuffmanBitCode>();
		//
		// int maxIndex = -1;
		// for (byte value : codes.keySet()) {
		// HuffmanBitCode code = codes.get(value);
		// code.value = value;
		//
		// if (code.bitLentgh == len) {
		// codesOfTheLength.add(code);
		// int index = (code.bitCode >>> prevLen) & (~(0 << len));
		// if (index > maxIndex)
		// maxIndex = index;
		// System.out.printf("\tfor code: code=%d\tlength=%d\tvalue=%d\n",
		// code.bitCode, code.bitLentgh,
		// code.value);
		// System.out.println("\tindex=" + index);
		// }
		// }
		//
		// node.leaf = new HuffmanBitCode[maxIndex + 1];
		//
		// for (HuffmanBitCode code : codesOfTheLength) {
		// int index = (code.bitCode >>> prevLen) & (~(0 << len));
		// node.leaf[index] = code;
		// System.out.printf("index=%d\tcode=%d\tlength=%d\tvalue=%d\n", index,
		// code.bitCode, code.bitLentgh,
		// code.value);
		// }
		//
		// prevLen = len;
		// }
	}

	private void printNode(Node node) {

	}

	private Node bn(int prevLen, Collection<HuffmanBitCode> codes) {
		Node node = new Node();

		int len = Integer.MAX_VALUE;
		for (HuffmanBitCode code : codes) {
			if (code.bitLentgh > prevLen && len > code.bitLentgh)
				len = code.bitLentgh;
		}

		node.len = len - prevLen;

		Map<Integer, List<HuffmanBitCode>> map = new HashMap<Integer, List<HuffmanBitCode>>();
		int maxLeafIndex = 0;
		int maxChildIndex = 0;

		for (HuffmanBitCode code : codes) {
			int index = (code.bitCode >>> (code.bitLentgh-len)) & (~(-1 << (len-prevLen)));

			if (code.bitLentgh == len)
				maxLeafIndex = maxLeafIndex > index ? maxLeafIndex : index;
			if (code.bitLentgh > len)
				maxChildIndex = maxChildIndex > index ? maxChildIndex : index;
		}

		node.leaf = new HuffmanBitCode[maxLeafIndex + 1];
		for (HuffmanBitCode code : codes) {
//			int index = code.bitCode & (~(-1 << (len-prevLen)));
			int index = (code.bitCode >>> (code.bitLentgh-len)) & (~(-1 << (len-prevLen)));

			if (code.bitLentgh == len) {
				node.leaf[index] = code;
//				System.out.println("Leaf");
			}
			if (code.bitLentgh > len) {
//				System.out.println("Node");
				List<HuffmanBitCode> list = map.get(index);
				if (list == null) {
					list = new ArrayList<HuffmanBitCode>();
					map.put(index, list);
				}
				list.add(code);
			}
//			System.out.printf("%d\t%d\t%d\n", code.bitCode, code.bitLentgh, code.value);
//			System.out.printf("\t%d\t%d\t%d\n", prevLen, len, index);
		}

		node.children = new Node[maxChildIndex + 1];
		for (Integer index : map.keySet()) {
			node.children[index] = bn(len, map.get(index));
		}

		return node;
	}

	@Override
	public Byte read(BitInputStream bis) throws IOException {
		Node node = root;
		while (node != null) {
			int readLongBits = bis.readBits(node.len);
			if (readLongBits < node.leaf.length) {
				HuffmanBitCode code = node.leaf[readLongBits];
				if (code != null)
					return code.value;
			}
			node = node.children[readLongBits];
		}
		throw new RuntimeException("Bit code not found. ");
	}

	@Override
	public long write(BitOutputStream bos, Byte symbol) throws IOException {
		HuffmanBitCode bitCode = bitCodes[symbol & 0xFF];
		if (bitCode == null)
			throw new RuntimeException("Huffman code not found for value: " + symbol);
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
		byte value;
	}

	private static <T> void getBitCode(HuffmanTree<Byte> tree, HuffmanBitCode code, Map<Byte, HuffmanBitCode> codes) {
		if (tree instanceof HuffmanLeaf) {
			HuffmanLeaf<Byte> leaf = (HuffmanLeaf<Byte>) tree;
			HuffmanBitCode readyCode = new HuffmanBitCode();
			readyCode.bitCode = code.bitCode;
			readyCode.bitLentgh = code.bitLentgh;
			readyCode.value = leaf.value;
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

	public void clear() {
		codes.clear();
		codes = null;
		tree = null;
	}

	public static void main(String[] args) throws IOException {
		System.out.println("Examples: ");
		printExample2("".getBytes());
		printExample2("A".getBytes());
		printExample2("AAAA".getBytes());
		printExample2("AB".getBytes());
		printExample2("AAB".getBytes());
		printExample2("AABB".getBytes());
		printExample2("AABBC".getBytes());
		printExample2("AAABB".getBytes());
		printExample2("AAABBC".getBytes());
		printExample2("AAABBBC".getBytes());
		printExample2("AAABBCD".getBytes());
		printExample2(("All the World's a Stage\n" + "\n" + "All the world's a stage,\n"
				+ "And all the men and women merely players;\n" + "They have their exits and their entrances,\n"
				+ "And one man in his time plays many parts,\n" + "His acts being seven ages. At first, the infant,\n"
				+ "Mewling and puking in the nurse's arms.\n" + "Then the whining schoolboy, with his satchel\n"
				+ "And shining morning face, creeping like snail\n" + "Unwillingly to school. And then the lover,\n"
				+ "Sighing like furnace, with a woeful ballad\n" + "Made to his mistress' eyebrow. Then a soldier,\n"
				+ "Full of strange oaths and bearded like the pard,\n"
				+ "Jealous in honor, sudden and quick in quarrel,\n" + "Seeking the bubble reputation\n"
				+ "Even in the cannon's mouth. And then the justice,\n"
				+ "In fair round belly with good capon lined,\n" + "With eyes severe and beard of formal cut,\n"
				+ "Full of wise saws and modern instances;\n" + "And so he plays his part. The sixth age shifts\n"
				+ "Into the lean and slippered pantaloon,\n" + "With spectacles on nose and pouch on side;\n"
				+ "His youthful hose, well saved, a world too wide\n"
				+ "For his shrunk shank, and his big manly voice,\n" + "Turning again toward childish treble, pipes\n"
				+ "And whistles in his sound. Last scene of all,\n" + "That ends this strange eventful history,\n"
				+ "Is second childishness and mere oblivion,\n"
				+ "Sans teeth, sans eyes, sans taste, sans everything. \n" + "William Shakespeare").getBytes());
	}

	private static void printExample2(byte[] data) throws IOException {
		ByteFrequencies bf = new ByteFrequencies();
		for (byte b : data)
			bf.add(b);

		int[] freqs = bf.getFrequencies();
		byte[] values = bf.getValues();

		HuffmanTree<Byte> tree2 = HuffmanCode.buildTree(freqs, Utils.autobox(values));

		System.out.println("For data: \"" + new String(data) + "\"");
		System.out.println("Data size in bits: " + data.length * 8);
		System.out.println("Values: " + Arrays.toString(values));
		System.out.println("Freqs: " + Arrays.toString(freqs));

		System.out.println("Tree: ");
		StringBuffer sb = new StringBuffer();
		HuffmanCode.printTree(tree2, sb, System.out);

		HuffmanByteCodec3 codec = new HuffmanByteCodec3(tree2);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);
		int len = 0;
		for (byte b : data)
			len += codec.write(bos, b);
		bos.close();

		System.out.println("Output size in bits: " + len);
		String bits = Utils.toBitString(baos.toByteArray());
		sb = new StringBuffer("\"");
		sb.append(bits.substring(0, len));
		sb.append(",");
		sb.append(bits.substring(len));
		sb.append("\"");
		System.out.println(sb.toString());
		System.out.println();

	}
}
