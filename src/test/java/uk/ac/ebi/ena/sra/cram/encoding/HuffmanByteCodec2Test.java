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

import static org.hamcrest.CoreMatchers.*;
import static org.junit.Assert.assertThat;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Random;

import org.junit.Ignore;
import org.junit.Test;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DebuggingBitOuputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

public class HuffmanByteCodec2Test {

	@Test
//	@Ignore
	public void test_write_1() throws IOException {
		Byte[] values = new Byte[] { 1, 2 };
		int[] charFreqs = new int[] { 1, 2 };
		HuffmanTree<Byte> tree = HuffmanCode.buildTree(charFreqs, values);

		HuffmanByteCodec2 codec = new HuffmanByteCodec2(tree);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		codec.write(bos, (byte) 1);
		bos.flush();
		byte[] buf = baos.toByteArray();

		assertThat(buf, equalTo(new byte[] { 0 }));
	}

	@Test
//	@Ignore
	public void test_write_2() throws IOException {
		Byte[] values = new Byte[] { 1, 2 };
		int[] charFreqs = new int[] { 1, 2 };
		HuffmanTree<Byte> tree = HuffmanCode.buildTree(charFreqs, values);

		HuffmanByteCodec2 codec = new HuffmanByteCodec2(tree);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		codec.write(bos, (byte) 2);
		bos.flush();
		byte[] buf = baos.toByteArray();

		assertThat(buf, equalTo(new byte[] { (byte) (1 << 7) }));
	}

	@Test
//	@Ignore
	public void test_write_1_2_3_4() throws IOException {
		Byte[] values = new Byte[] { 1, 2, 3, 4 };
		int[] charFreqs = new int[] { 1, 2, 3, 4 };
		HuffmanTree<Byte> tree = HuffmanCode.buildTree(charFreqs, values);

		HuffmanByteCodec2 codec = new HuffmanByteCodec2(tree);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		// 110:
		codec.write(bos, (byte) 1);

		// 111:
		codec.write(bos, (byte) 2);

		// 1:
		codec.write(bos, (byte) 3);

		// 0:
		codec.write(bos, (byte) 4);

		bos.flush();
		byte[] buf = baos.toByteArray();

		// 1101111000000000:
		assertThat(buf, equalTo(new byte[] { (byte) (1 << 7 | 1 << 6 | 1 << 4 | 1 << 3 | 1 << 2 | 1 << 1), 0 }));
	}

	@Test
//	@Ignore
	public void test_write_read_random() throws IOException {
		int maxTests = 1000;
		Random random = new Random();
		Byte[] alphabet = new Byte[256];
		int[] charFreqs = new int[alphabet.length];
		for (int b = Byte.MIN_VALUE; b <= Byte.MAX_VALUE; b++) {
			alphabet[b & 0xFF] = (byte) b;
			charFreqs[b & 0xFF] = b + 128;
		}

		Byte[] values = new Byte[maxTests];
		for (int i = 0; i < maxTests; i++) {
			values[i] = (byte) alphabet.length;
			for (int j = 0; j < random.nextInt(alphabet.length); j++)
				if (random.nextBoolean())
					values[i]--;
		}

		HuffmanTree<Byte> tree = HuffmanCode.buildTree(charFreqs, alphabet);

		HuffmanByteCodec2 codec = new HuffmanByteCodec2(tree);
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

	@Test//(timeout = 650)
	public void benchmark_write_read_random() throws IOException {
		int maxTests = 1000000;
		Random random = new Random();
		Byte[] alphabet = new Byte[] { 1, 2, 3, 4 };
		Byte[] values = new Byte[maxTests];
		for (int i = 0; i < maxTests; i++) {
			values[i] = (byte) alphabet.length;
			for (int j = 0; j < random.nextInt(alphabet.length); j++)
				if (random.nextBoolean())
					values[i]--;
		}

		int[] charFreqs = new int[] { 1, 2, 3, 4 };
		HuffmanTree<Byte> tree = HuffmanCode.buildTree(charFreqs, alphabet);
		BitCodec<Byte> codec = new HuffmanByteCodec2(tree);

		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);
//		BitOutputStream bos = new DebuggingBitOuputStream(System.out, '\n', new DefaultBitOutputStream(baos));

		for (int i = 0; i < maxTests; i++)
			codec.write(bos, values[i]);

		bos.flush();
		byte[] buf = baos.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);

		for (int i = 0; i < maxTests; i++)
			assertThat(codec.read(bis), equalTo(values[i]));
	}

	@Test
	public void benchmark_read_random() throws IOException {
		int maxTests = 10000;
		Byte[] alphabet = new Byte[] { 1, 2, 3, 4 };
		int[] charFreqs = new int[] { 1, 2, 3, 4 };
		HuffmanTree<Byte> tree = HuffmanCode.buildTree(charFreqs, alphabet);
		HuffmanByteCodec2 codec = new HuffmanByteCodec2(tree);
		byte[] buf = "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"
				.getBytes();

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);

		for (int i = 0; i < maxTests; i++) {
			Byte b = codec.read(bis);
			assertThat(b >= 1 && b <= 4, is(true));
		}

	}
}
