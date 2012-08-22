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
import java.util.Arrays;
import java.util.Random;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ExponentialDistributionImpl;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.junit.Ignore;
import org.junit.Test;

import uk.ac.ebi.ena.sra.compression.ArrayCompressionReport;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

public class GolombRiceCodecTest {

	@Test
	public void test_numberOfBits() {
		GolombRiceCodec codec = new GolombRiceCodec(2);
		for (long value = 0; value < 1000; value++)
			assertThat(codec.numberOfBits(value), is(value / 4 + 3));
	}

	@Test
	@Ignore
	public void printCodes_1_to_256() throws IOException {
		int golombRiceLogM = 2;
		GolombRiceCodec codec = new GolombRiceCodec(golombRiceLogM);
		for (int i = 0; i < 256; i++) {
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
			BitOutputStream bos = new DefaultBitOutputStream(baos);
			int len = (int) codec.write(bos, (long) i);
			bos.flush();

			byte[] buf = baos.toByteArray();

			ByteArrayInputStream bais = new ByteArrayInputStream(buf);
			BitInputStream bis = new DefaultBitInputStream(bais);
			long number = codec.read(bis);
			System.out.printf("%d: %d\t%s\t%d\t%s\n", i, number, Utils.toBitString(buf).subSequence(0, len), len,
					Utils.toBitString(buf));
		}
	}

	@Test
	public void test_0() throws IOException {
		long value = 0;
		long bitsLen = 3;
		int golombRiceLogM = 2;
		GolombRiceCodec codec = new GolombRiceCodec(golombRiceLogM);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);
		long len = codec.write(bos, value);
		bos.flush();

		assertThat(len, is(bitsLen));

		byte[] buf = baos.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);
		long number = codec.read(bis);
		assertThat(number, is(value));

		String bitsString = Utils.toBitString(buf);
		assertThat(bitsString, equalTo("10000000"));

		String cutBitsString = bitsString.substring(0, (int) len);
		assertThat(cutBitsString, equalTo("100"));
	}

	@Test
	public void test_1() throws IOException {
		long value = 1;
		long bitsLen = 3;
		int golombRiceLogM = 2;
		GolombRiceCodec codec = new GolombRiceCodec(golombRiceLogM);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);
		long len = codec.write(bos, value);
		bos.flush();

		assertThat(len, is(bitsLen));

		byte[] buf = baos.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);
		long number = codec.read(bis);
		assertThat(number, is(value));

		String bitsString = Utils.toBitString(buf);
		assertThat(bitsString, equalTo("10100000"));

		String cutBitsString = bitsString.substring(0, (int) len);
		assertThat(cutBitsString, equalTo("101"));
	}

	@Test
	public void test_20() throws IOException {
		long value = 20;
		long bitsLen = 8;
		int golombRiceLogM = 2;
		GolombRiceCodec codec = new GolombRiceCodec(golombRiceLogM);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);
		long len = codec.write(bos, value);
		bos.flush();

		assertThat(len, is(bitsLen));

		byte[] buf = baos.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);
		long number = codec.read(bis);
		assertThat(number, is(value));

		String bitsString = Utils.toBitString(buf);
		assertThat(bitsString, equalTo("00000100"));

		String cutBitsString = bitsString.substring(0, (int) len);
		assertThat(cutBitsString, equalTo("00000100"));
	}

	@Test
	public void test_255() throws IOException {
		long value = 255;
		long bitsLen = 66;
		int golombRiceLogM = 2;
		GolombRiceCodec codec = new GolombRiceCodec(golombRiceLogM);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);
		long len = codec.write(bos, value);
		bos.flush();

		assertThat(len, is(bitsLen));

		byte[] buf = baos.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);
		long number = codec.read(bis);
		assertThat(number, is(value));

		String bitsString = Utils.toBitString(buf);
		assertThat(bitsString, equalTo("000000000000000000000000000000000000000000000000000000000000000111000000"));

		String cutBitsString = bitsString.substring(0, (int) len);
		assertThat(cutBitsString, equalTo("000000000000000000000000000000000000000000000000000000000000000111"));
	}

	@Test(timeout = 3000)
	public void benchmark_uniform_distr() throws IOException {
		int maxNumbers = 1000000;
		for (int log2m = 1; log2m < 6; log2m++) {
			int golombParameter = 1 << log2m;
			GolombRiceCodec codec = new GolombRiceCodec(log2m);
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
			BitOutputStream bos = new DefaultBitOutputStream(baos);

			Random random = new Random();
			int maxRandom = golombParameter * 2;
			for (int i = 0; i < maxNumbers; i++)
				codec.write(bos, new Long(random.nextInt(maxRandom)));

			bos.flush();
			byte[] buf = baos.toByteArray();
			assertThat(buf.length > 300000, is(true));
			assertThat(buf.length < 900000, is(true));

			ArrayCompressionReport report = new ArrayCompressionReport("compressability");
			report.run(buf);

			float ratioToBzip = (float) buf.length / report.getBzip2Size();
			assertThat(ratioToBzip < 1.5, is(true));
			assertThat(ratioToBzip > 1F, is(true));

			float ratioToGzip = (float) buf.length / report.getGzipSize();
			assertThat(ratioToGzip < 1.5, is(true));
			assertThat(ratioToGzip > 1F, is(true));

			// System.out
			// .printf("Universal distribution (average=%d), golomb param=%d, bits per number: %2.2f, %s\n",
			// maxRandom / 2, golombParameter, 8f * baos.size()
			// / maxNumbers, report.toString());

			baos.close();
		}
	}

	@Test(timeout = 3600)
	public void benchmark_exponential_distr() throws IOException, MathException {
		int maxNumbers = 1000000;
		for (int log2m = 1; log2m < 6; log2m++) {
			int golombParameter = 1 << log2m;
			GolombRiceCodec codec = new GolombRiceCodec(log2m);
			ExponentialDistributionImpl expDistr = new ExponentialDistributionImpl(golombParameter);
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
			BitOutputStream bos = new DefaultBitOutputStream(baos);

			for (int i = 0; i < maxNumbers; i++)
				codec.write(bos, (long) expDistr.sample());

			bos.flush();
			byte[] buf = baos.toByteArray();
			assertThat(buf.length > 300000, is(true));
			assertThat(buf.length < 900000, is(true));

			ArrayCompressionReport report = new ArrayCompressionReport("compressability");
			report.run(buf);

			float ratioToBzip = (float) buf.length / report.getBzip2Size();
			assertThat(ratioToBzip < 1.11, is(true));
			assertThat(ratioToBzip > 0.99F, is(true));

			float ratioToGzip = (float) buf.length / report.getGzipSize();
			assertThat(ratioToGzip < 1.11, is(true));
			assertThat(ratioToGzip > 0.99F, is(true));

			baos.close();
		}
	}

	@Test(timeout = 11000)
	public void benchmark_normal_distr() throws IOException, MathException {
		int maxNumbers = 1000000;
		for (int log2m = 1; log2m < 6; log2m++) {
			int golombParameter = 1 << log2m;
			GolombRiceCodec codec = new GolombRiceCodec(log2m);
			double mean = golombParameter * 10;
			double sd = golombParameter * 5;
			NormalDistributionImpl expDistr = new NormalDistributionImpl(mean, sd);
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
			BitOutputStream bos = new DefaultBitOutputStream(baos);

			for (int i = 0; i < maxNumbers; i++) {
				double sample = 0;
				while ((sample = expDistr.sample()) < 0)
					;
				codec.write(bos, (long) Math.abs(expDistr.sample()));
			}

			bos.flush();
			byte[] buf = baos.toByteArray();
			assertThat(buf.length > 300000, is(true));
			assertThat(buf.length < 2000000, is(true));

			ArrayCompressionReport report = new ArrayCompressionReport("compressability");
			report.run(baos.toByteArray());

			// System.out
			// .printf("Normal distribution (average=%.2f, sd=%.2f), golomb param=%d, bits per number: %2.2f, %s\n",
			// mean, sd, golombParameter, 8f * baos.size()
			// / maxNumbers, report.toString());

			float ratioToBzip = (float) buf.length / report.getBzip2Size();
			assertThat(ratioToBzip < 2, is(true));
			assertThat(ratioToBzip > 1F, is(true));

			float ratioToGzip = (float) buf.length / report.getGzipSize();
			assertThat(ratioToGzip < 2, is(true));
			assertThat(ratioToGzip > 1F, is(true));

			baos.close();
		}
	}

	@Test(timeout = 300)
	public void becnmark_Write() throws IOException {
		int maxNumbers = 1000000;
		int log2m = 2;
		GolombRiceCodec codec = new GolombRiceCodec(log2m);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		for (int i = 0; i < maxNumbers; i++)
			codec.write(bos, 20L);

		bos.flush();
		baos.close();
	}

	@Test(timeout = 200)
	public void becnmark_Read() throws IOException {
		int maxNumbers = 1000000;
		int log2m = 2;
		GolombRiceCodec codec = new GolombRiceCodec(log2m);
		byte oneByteValue = 20;

		byte[] buf = new byte[maxNumbers];
		Arrays.fill(buf, oneByteValue);
		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);
		for (int i = 0; i < maxNumbers; i++)
			codec.read(bis);
	}
}
