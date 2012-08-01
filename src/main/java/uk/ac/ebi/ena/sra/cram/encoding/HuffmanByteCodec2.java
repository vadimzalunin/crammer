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
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
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

public class HuffmanByteCodec2 implements BitCodec<Byte> {
	private HuffmanTree<Byte> tree;
	private TreeMap<Byte, HuffmanBitCode> codes;
	private HuffmanBitCode[] bitCodes = new HuffmanBitCode[256];
	// assuming code length cannot be more than 1024:
	private Integer[] codeLentghSorted;
	private Map<Integer, Map<Long, Byte>> codeCache = new HashMap<Integer, Map<Long, Byte>>();
	private Map<Long, Byte>[] codeMaps ;
	
	private static long readCounter = 0 ;
	private static long readTime = 0 ;
	
	private static long readLongCounter = 0 ;
	private static long readLongTime = 0 ;
	
	public static void dump (){
//		System.out.println("HuffmanByteCodec2 read stats: count=" + readCounter + "; millis=" + readTime);
//		System.out.println("HuffmanByteCodec2 read long bits stats: count=" + readLongCounter + "; millis=" + readLongTime);
	}

	public static HuffmanByteCodec2 build(ByteFrequencies bf) {
		HuffmanTree<Byte> tree = HuffmanCode.buildTree(bf.getFrequencies(), Utils.autobox(bf.getValues()));
		return new HuffmanByteCodec2(tree);
	}
	
	public HuffmanByteCodec2(HuffmanTree<Byte> tree) {
		super();
		this.tree = tree;
		codes = new TreeMap<Byte, HuffmanBitCode>();
		getBitCode(tree, new HuffmanBitCode(), codes);
		
		if (codes.isEmpty()) return ;

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
		
		codeMaps = new Map[codeLentghSorted[codeLentghSorted.length-1]+1] ;
		for (int len:codeLentghSorted) {
			codeMaps[len] = codeCache.get(len) ;
		}
	}

	@Override
	public Byte read(BitInputStream bis) throws IOException {
//		long time = System.currentTimeMillis() ;
		long buf = 0;
		int bitsRead = 0;
		for (int len : codeLentghSorted) {
			buf = buf << (len - bitsRead);
			
//			readLongCounter ++ ;
//			long rlTime = System.currentTimeMillis() ;
			long readLongBits = bis.readLongBits(len - bitsRead) ;
//			readLongTime += System.currentTimeMillis() -rlTime ;
			
			buf = buf | readLongBits;
			
			
			bitsRead = len;
			Map<Long, Byte> codeMap = codeMaps[len];
//			if (codeMap == null)
//				continue;
			Byte result = codeMap.get(buf);
			if (result != null) {
//				readTime += System.currentTimeMillis() - time ;
//				readCounter++ ;
				return result ;
			}
		}
		throw new RuntimeException("Bit code not found. Current state: " + bitsRead + " bits read, buf=" + buf);
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
}
