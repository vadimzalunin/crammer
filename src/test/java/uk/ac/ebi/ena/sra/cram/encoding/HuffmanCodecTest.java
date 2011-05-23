package uk.ac.ebi.ena.sra.cram.encoding;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.junit.Assert.assertThat;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Random;

import org.junit.Test;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

public class HuffmanCodecTest {

	@Test
	public void test_write_1() throws IOException {
		Integer[] values = new Integer[] { 1, 2 };
		int[] charFreqs = new int[] { 1, 2 };
		HuffmanTree<Integer> tree = HuffmanCode.buildTree(charFreqs, values);

		HuffmanCodec<Integer> codec = new HuffmanCodec<Integer>(tree);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		codec.write(bos, 1);
		bos.flush();
		byte[] buf = baos.toByteArray();

		assertThat(buf, equalTo(new byte[] { 0 }));
	}

	@Test
	public void test_write_2() throws IOException {
		Integer[] values = new Integer[] { 1, 2 };
		int[] charFreqs = new int[] { 1, 2 };
		HuffmanTree<Integer> tree = HuffmanCode.buildTree(charFreqs, values);

		HuffmanCodec<Integer> codec = new HuffmanCodec<Integer>(tree);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		codec.write(bos, 2);
		bos.flush();
		byte[] buf = baos.toByteArray();

		assertThat(buf, equalTo(new byte[] { (byte) (1 << 7) }));
	}

	@Test
	public void test_write_1_2_3_4() throws IOException {
		Integer[] values = new Integer[] { 1, 2, 3, 4 };
		int[] charFreqs = new int[] { 1, 2, 3, 4 };
		HuffmanTree<Integer> tree = HuffmanCode.buildTree(charFreqs, values);

		HuffmanCodec<Integer> codec = new HuffmanCodec<Integer>(tree);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		// 110:
		codec.write(bos, 1);

		// 111:
		codec.write(bos, 2);

		// 1:
		codec.write(bos, 3);

		// 0:
		codec.write(bos, 4);

		bos.flush();
		byte[] buf = baos.toByteArray();

		// 1101111000000000:
		assertThat(
				buf,
				equalTo(new byte[] {
						(byte) (1 << 7 | 1 << 6 | 1 << 4 | 1 << 3 | 1 << 2 | 1 << 1),
						0 }));
	}

	@Test
	public void test_write_read_random() throws IOException {
		int maxTests = 1000;
		Random random = new Random();
		Integer[] alphabet = new Integer[] { 1, 2, 3, 4 };
		Integer[] values = new Integer[maxTests];
		for (int i = 0; i < maxTests; i++) {
			values[i] = alphabet.length;
			for (int j = 0; j < random.nextInt(alphabet.length); j++)
				if (random.nextBoolean())
					values[i]--;
		}

		int[] charFreqs = new int[] { 1, 2, 3, 4 };
		HuffmanTree<Integer> tree = HuffmanCode.buildTree(charFreqs, alphabet);

		HuffmanCodec<Integer> codec = new HuffmanCodec<Integer>(tree);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		for (int i = 0; i < maxTests; i++)
			codec.write(bos, values[i]);

		bos.flush();
		byte[] buf = baos.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);

		for (int i = 0; i < maxTests; i++)
			assertThat(codec.read(bis), equalTo(values[i]));

	}
	
	@Test (timeout=450)
	public void benchmark_write_read_random() throws IOException {
		int maxTests = 1000000;
		Random random = new Random();
		Integer[] alphabet = new Integer[] { 1, 2, 3, 4 };
		Integer[] values = new Integer[maxTests];
		for (int i = 0; i < maxTests; i++) {
			values[i] = alphabet.length;
			for (int j = 0; j < random.nextInt(alphabet.length); j++)
				if (random.nextBoolean())
					values[i]--;
		}

		int[] charFreqs = new int[] { 1, 2, 3, 4 };
		HuffmanTree<Integer> tree = HuffmanCode.buildTree(charFreqs, alphabet);

		HuffmanCodec<Integer> codec = new HuffmanCodec<Integer>(tree);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		for (int i = 0; i < maxTests; i++)
			codec.write(bos, values[i]);

		bos.flush();
		byte[] buf = baos.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);

		for (int i = 0; i < maxTests; i++)
			assertThat(codec.read(bis), equalTo(values[i]));

	}
}
