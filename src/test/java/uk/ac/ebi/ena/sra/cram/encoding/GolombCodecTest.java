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

import org.junit.Ignore;
import org.junit.Test;

import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

public class GolombCodecTest {

	@Test
	public void test_numberOfBits() {
		GolombCodec codec = new GolombCodec(3);
		assertThat(codec.numberOfBits(1L), is(3L));
		assertThat(codec.numberOfBits(4L), is(4L));
		assertThat(codec.numberOfBits(8L), is(5L));
		assertThat(codec.numberOfBits(10L), is(6L));
	}

	@Test
	public void test_write_14_3() throws IOException {
		int m = 3;
		GolombCodec codec = new GolombCodec(m);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);
		long value = 14L;
		int len = (int) codec.write(bos, value);
		assertThat(len, is(7));
		bos.flush();

		byte[] buf = baos.toByteArray();
		assertThat(Utils.toBitString(buf).substring(0, len), equalTo("1111011"));
	}

	@Test
	public void test_read_14_3() throws IOException {
		int m = 3;
		long value = 14L;
		// 1111011:
		byte[] buf = new byte[] { (byte) (123 << 1) };

		GolombCodec codec = new GolombCodec(m);
		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);

		Long readValue = codec.read(bis);
		assertThat(readValue, is(value));
	}

	@Test
	public void test_write_14_4() throws IOException {
		int m = 4;
		GolombCodec codec = new GolombCodec(m);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);
		long value = 14L;
		int len = (int) codec.write(bos, value);
		assertThat(len, is(6));
		bos.flush();

		byte[] buf = baos.toByteArray();
		assertThat(Utils.toBitString(buf).substring(0, len), equalTo("111010"));
	}

	@Test
	public void test_read_14_43() throws IOException {
		int m = 4;
		long value = 14L;
		// 111010:
		byte[] buf = new byte[] { (byte) (58 << 2) };

		GolombCodec codec = new GolombCodec(m);
		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);

		Long readValue = codec.read(bis);
		assertThat(readValue, is(value));
	}

	@Test(expected = IllegalArgumentException.class)
	public void test_write_read_1() throws IOException {
		new GolombCodec(1);
	}

	@Test
	@Ignore
	public void printCodes_1_to_256() throws IOException {
		int m = 3;
		GolombCodec codec = new GolombCodec(m);
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
}
