package uk.ac.ebi.ena.sra.compression;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ExponentialDistributionImpl;
import org.apache.commons.math.stat.Frequency;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanLeaf;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanNode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

import com.colloquial.arithcode.AdaptiveUnigramModel;
import com.colloquial.arithcode.ArithCodeOutputStream;
import com.colloquial.arithcode.ArithEncoder;
import com.colloquial.arithcode.BitOutput;

public class TestArithmeticCode {

	public static void main(String[] args) throws IOException, MathException {
		AdaptiveUnigramModel modelIn = new AdaptiveUnigramModel();
		// PPMModel modelIn = new PPMModel(1);
		// UniformModel modelIn = UniformModel.MODEL ;
		ByteArrayOutputStream arithCodeBaos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(arithCodeBaos) ;
		ArithEncoder encoder = new ArithEncoder(bos) ;
		ArithCodeOutputStream arithCodeOutputStream = new ArithCodeOutputStream(encoder, modelIn);

		ByteArrayOutputStream gzipBaos = new ByteArrayOutputStream();
		GZIPOutputStream gzipOS = new GZIPOutputStream(gzipBaos);

		long testSize = 1000000L;

		ExponentialDistributionImpl expDistr = new ExponentialDistributionImpl(
				0.001);
		Frequency freq = new Frequency();
		Random random = new Random() ;
		for (long i = 0; i < testSize; i++) {
//			byte randomValue = (byte) (255 * expDistr.sample());
			byte randomValue = 0 ;
			switch (random.nextInt(100)) {
			case 0:
				randomValue = 4 ;
				break;

			default:
				randomValue = (byte) random.nextInt(4) ;
				break;
			}
			freq.addValue(randomValue);
			byte[] buf = new byte[] { randomValue };
			arithCodeOutputStream.write(buf);
			gzipOS.write(buf);
		}
		arithCodeOutputStream.close();
		gzipOS.close();

		System.out.println("Arithmetic code: " + arithCodeBaos.size());
		System.out.println("Gzip code: " + gzipBaos.size());
		System.out.println("Bases: " + testSize);

		System.out.println(freq.toString());

		byte[] alphabet = new byte[freq.getUniqueCount()];
		Iterator<Comparable<?>> valuesIterator = freq.valuesIterator();
		int i = 0;
		int[] frequencies = new int[33 + freq.getUniqueCount()];
		while (valuesIterator.hasNext()) {
			long b = (Long) valuesIterator.next();
			alphabet[i] = (byte) b;
			frequencies[i + 33] = (int) freq.getCount(b);
			i++;
		}

		Long[] values = new Long[frequencies.length] ;
		for (int v=0; v<values.length; v++) values[v] = (long) v ;
		HuffmanTree<Long> huffmanTree = HuffmanCode.buildTree(frequencies, values);
		System.out.println();

		ArrayList<HuffmanBitCode> list = new ArrayList<HuffmanBitCode>();
		getBitCode(huffmanTree, new HuffmanBitCode(), list);
		Collections.sort(list, new Comparator<HuffmanBitCode>() {

			@Override
			public int compare(HuffmanBitCode o1, HuffmanBitCode o2) {
				return o2.frequency - o1.frequency;
			}
		});

		long length = 0L;
		for (HuffmanBitCode code : list) {
			String binaryCode = String.format("%1$#" + code.bitLentgh + "s",
					Long.toBinaryString(code.bitCode), code.bitLentgh);
			System.out.println((char) code.value + "\t" + code.frequency + "\t"
					+ binaryCode.replaceAll(" ", "0"));
			length += code.frequency * code.bitLentgh;
		}
		System.out.println(1 + length / 8);

	}

	private static class HuffmanBitCode {
		private long bitCode;
		private int bitLentgh;
		private long value;
		private int frequency;
	}

	private static void getBitCode(HuffmanTree<Long> tree, HuffmanBitCode code,
			List<HuffmanBitCode> codes) {
		assert tree != null;
		if (tree instanceof HuffmanLeaf) {
			HuffmanLeaf<Long> leaf = (HuffmanLeaf<Long>) tree;
			HuffmanBitCode readyCode = new HuffmanBitCode();
			readyCode.bitCode = code.bitCode;
			readyCode.bitLentgh = code.bitLentgh;
			readyCode.frequency = leaf.frequency;
			readyCode.value = leaf.value;
			codes.add(readyCode);
			return;

		} else if (tree instanceof HuffmanNode) {
			HuffmanNode<Long> node = (HuffmanNode<Long>) tree;

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
}
