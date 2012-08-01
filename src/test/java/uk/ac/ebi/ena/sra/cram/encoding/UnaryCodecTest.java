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

import org.junit.Test;

import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

public class UnaryCodecTest {

	@Test
	public void test_numberOfBits() {
		UnaryCodec codec = new UnaryCodec(true, 0);
		for (long value = 0; value < 1000; value++)
			assertThat(codec.numberOfBits(value), is(value + 1));
	}

	@Test
	public void test_write_read_0() throws IOException {
		UnaryCodec codec = new UnaryCodec(false, 0);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		long value = 0L;
		long len = codec.write(bos, value);
		assertThat(len, is(value + 1));

		bos.flush();
		byte[] buf = baos.toByteArray();
		String bitString = Utils.toBitString(buf);

		assertThat(bitString, equalTo("00000000"));

		BitInputStream bis = new DefaultBitInputStream(
				new ByteArrayInputStream(buf));
		Long readValue = codec.read(bis);
		assertThat(readValue, equalTo(value));
	}

	@Test
	public void test_write_read_20() throws IOException {
		UnaryCodec codec = new UnaryCodec(false, 0);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		long value = 20L;
		long len = codec.write(bos, value);
		assertThat(len, is(value + 1));

		bos.flush();
		byte[] buf = baos.toByteArray();
		String bitString = Utils.toBitString(buf);

		assertThat(bitString, equalTo("111111111111111111110000"));

		BitInputStream bis = new DefaultBitInputStream(
				new ByteArrayInputStream(buf));
		Long readValue = codec.read(bis);
		assertThat(readValue, equalTo(value));
	}

	@Test
	public void test_write_read_0_to_1000() throws IOException {
		UnaryCodec codec = new UnaryCodec(false, 0);

		for (long i = 0; i < 10000; i++) {
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
			BitOutputStream bos = new DefaultBitOutputStream(baos);
			long len = codec.write(bos, i);
			assertThat(len, is(i + 1));

			bos.flush();
			byte[] buf = baos.toByteArray();
			assertThat("i=" + i + "; " + "buf len=" + buf.length + "; len="
					+ len, (long) buf.length, is(1 + (len - 1 >> 3)));

			BitInputStream bis = new DefaultBitInputStream(
					new ByteArrayInputStream(buf));
			Long readValue = codec.read(bis);
			assertThat(readValue, equalTo(i));
		}
	}

	@Test(timeout = 1000)
	public void benchmark_write_read() throws IOException {
		long maxValues = 1000000;

		UnaryCodec codec = new UnaryCodec(false, 0);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		for (long i = 0; i < maxValues; i++)
			codec.write(bos, 20L);

		bos.flush();
		byte[] buf = baos.toByteArray();

		BitInputStream bis = new DefaultBitInputStream(
				new ByteArrayInputStream(buf));

		for (long i = 0; i < maxValues; i++) {
			codec.read(bis);
		}
	}
}
