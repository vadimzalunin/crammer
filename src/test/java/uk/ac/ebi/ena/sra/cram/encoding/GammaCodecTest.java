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

import org.apache.commons.math.util.MathUtils;
import org.junit.Test;

import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

public class GammaCodecTest {
	@Test
	public void test_numberOfBits() {
		GammaCodec codec = new GammaCodec(0);

		assertThat(codec.numberOfBits((long) 1), is(1L));
		assertThat(codec.numberOfBits((long) 2), is(3L));

		for (long value = 2; value < 1000; value++) {
			int betaCodeLength = 1 + (int) MathUtils.log(2, value);
			assertThat(codec.numberOfBits(value), is((long) betaCodeLength * 2 - 1));
		}
	}

	@Test
	public void test_write_read_1() throws IOException {
		GammaCodec codec = new GammaCodec(0);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		long value = 1L;
		codec.write(bos, value);

		bos.flush();
		byte[] buf = baos.toByteArray();
		String bitString = Utils.toBitString(buf);

		assertThat(bitString, equalTo("10000000"));

		BitInputStream bis = new DefaultBitInputStream(new ByteArrayInputStream(buf));
		Long readValue = codec.read(bis);
		assertThat(readValue, equalTo(value));
	}

	@Test
	public void test_write_read_2() throws IOException {
		GammaCodec codec = new GammaCodec(0);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		long value = 2L;
		codec.write(bos, value);

		bos.flush();
		byte[] buf = baos.toByteArray();
		String bitString = Utils.toBitString(buf);

		assertThat(bitString, equalTo("01000000"));

		BitInputStream bis = new DefaultBitInputStream(new ByteArrayInputStream(buf));
		Long readValue = codec.read(bis);
		assertThat(readValue, equalTo(value));
	}

	@Test
	public void test_write_read_20() throws IOException {
		GammaCodec codec = new GammaCodec(0);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		long value = 20L;
		codec.write(bos, value);

		bos.flush();
		byte[] buf = baos.toByteArray();
		String bitString = Utils.toBitString(buf);

		// 20 is 10100 in unary, prepend it with 4 zeros and append with zero's
		// to pad to whole bytes:
		assertThat(bitString, equalTo("0000101000000000"));

		BitInputStream bis = new DefaultBitInputStream(new ByteArrayInputStream(buf));
		Long readValue = codec.read(bis);
		assertThat(readValue, equalTo(value));
	}

	@Test
	public void test_write_read_128() throws IOException {
		GammaCodec codec = new GammaCodec(0);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		long value = 128L;
		codec.write(bos, value);

		bos.flush();
		byte[] buf = baos.toByteArray();
		String bitString = Utils.toBitString(buf);

		// 128 is 10000000 in unary, prepend it with 7 zeros and append with
		// zero's
		// to pad to whole bytes:
		assertThat(bitString, equalTo("0000000100000000"));

		BitInputStream bis = new DefaultBitInputStream(new ByteArrayInputStream(buf));
		Long readValue = codec.read(bis);
		assertThat(readValue, equalTo(value));
	}

	@Test
	public void test_write_read_1_to_1000_offset_0() throws IOException {
		GammaCodec codec = new GammaCodec(0);

		for (long i = 1; i < 1000; i++) {
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
			BitOutputStream bos = new DefaultBitOutputStream(baos);
			codec.write(bos, i);

			bos.flush();
			byte[] buf = baos.toByteArray();

			BitInputStream bis = new DefaultBitInputStream(new ByteArrayInputStream(buf));
			Long readValue = codec.read(bis);
			assertThat(readValue, equalTo(i));
		}
	}

	@Test
	public void test_write_read_0_to_1000_offset_1() throws IOException {
		GammaCodec codec = new GammaCodec(1);

		for (long i = 0; i < 1000; i++) {
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
			BitOutputStream bos = new DefaultBitOutputStream(baos);
			codec.write(bos, i);

			bos.flush();
			byte[] buf = baos.toByteArray();

			BitInputStream bis = new DefaultBitInputStream(new ByteArrayInputStream(buf));
			Long readValue = codec.read(bis);
			assertThat(readValue, equalTo(i));
		}
	}

	@Test(timeout = 1000)
	public void benchmark_write_read() throws IOException {
		long maxValues = 1000000;

		GammaCodec codec = new GammaCodec(1);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		for (long i = 0; i < maxValues; i++)
			codec.write(bos, 20L);

		bos.flush();
		byte[] buf = baos.toByteArray();

		BitInputStream bis = new DefaultBitInputStream(new ByteArrayInputStream(buf));

		for (long i = 0; i < maxValues; i++) {
			codec.read(bis);
		}
	}

}
