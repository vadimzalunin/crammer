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
package uk.ac.ebi.ena.sra.cram.io;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.Arrays;

import org.junit.Test;

import uk.ac.ebi.ena.sra.cram.Utils;

public class BitInputStreamTest {

	@Test
	public void test_ReadBits_8_bits() throws IOException {
		byte value = 4;
		byte[] buf = new byte[] { value };
		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);
		int readBits = bis.readBits(8);

		assertThat((byte) readBits, is(value));
	}

	@Test
	public void test_ReadBits_3_bits_of_00000100() throws IOException {
		byte value = 4;
		byte[] buf = new byte[] { value };
		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);
		int readBits = bis.readBits(3);

		assertThat(readBits, is(0));
	}

	@Test
	public void test_ReadBits_3_bits_of_01000000() throws IOException {
		byte value = (byte) (1 << 7);
		byte[] buf = new byte[] { value };
		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);
		int readBits = bis.readBits(3);

		assertThat(readBits, is(4));
	}

	@Test
	public void test_ReadBits_31_bits_of_00000000000000000000000000000010()
			throws IOException {
		// 00000000000000000000000000000010:
		int value = 1 << 1;
		byte[] buf = Utils.toBytes(value);
		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);

		// 00000000000000000000000000000001:
		int readBits = bis.readBits(31);

		assertThat(readBits, is(1));
	}

	@Test
	public void test_ReadBits_int_by_bit() throws IOException {
		int setBit = 2;
		int value = 1 << setBit;
		byte[] buf = Utils.toBytes(value);
		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);
		for (int bit = 31; bit > -1; bit--) {
			int readBits = bis.readBits(1);

			if (bit == setBit)
				assertThat(readBits, is(1));
			else
				assertThat(readBits, is(0));
		}
	}

	@Test(timeout = 350)
	public void becnhmar_ReadBits_32() throws IOException {
		int maxValues = 1000000;
		byte[] buf = new byte[maxValues * 4];
		byte fillByte = (byte) (64 + 16 + 4 + 1);
		Arrays.fill(buf, fillByte);
		int value = 0;
		for (int i = 0; i < 4; i++)
			value = (value << 8) + fillByte;
		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);
		for (int i = 0; i < maxValues; i++)
			bis.readBits(32);

	}

	@Test(timeout = 550)
	public void becnhmar_ReadLongBits_64() throws IOException {
		int maxValues = 1000000;
		byte[] buf = new byte[maxValues * 8];
		byte fillByte = (byte) (64 + 16 + 4 + 1);
		Arrays.fill(buf, fillByte);
		int value = 0;
		for (int i = 0; i < 8; i++)
			value = (value << 8) + fillByte;
		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);
		for (int i = 0; i < maxValues; i++)
			bis.readLongBits(64);
	}

	@Test(timeout = 750)
	public void becnhmar_ReadLongBits_32() throws IOException {
		int maxValues = 1000000;
		byte[] buf = new byte[maxValues * 4];
		byte fillByte = (byte) (64 + 16 + 4 + 1);
		Arrays.fill(buf, fillByte);
		int value = 0;
		for (int i = 0; i < 4; i++)
			value = (value << 8) + fillByte;
		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		DefaultBitInputStream bis = new DefaultBitInputStream(bais);
		for (int i = 0; i < maxValues; i++)
			bis.readLongBits(32);
	}

	@Test
	public void test_ReadBits_long_3_bits_of_01000000() throws IOException {
		byte value = (byte) (1 << 7);
		byte[] buf = new byte[] { value };
		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		DefaultBitInputStream bis = new DefaultBitInputStream(bais);
		long readBits = bis.readLongBits(3);

		assertThat(readBits, is(4L));
	}

	@Test
	public void test_ReadBits_int_range() throws IOException {
		for (int i = 1; i < 1000000; i++) {
			int len = (int) (1 + Math.log(i) / Math.log(2));
			byte[] buf = Utils.toBytes(i << 32 - len);
			ByteArrayInputStream bais = new ByteArrayInputStream(buf);
			DefaultBitInputStream bis = new DefaultBitInputStream(bais);
			int readBits = bis.readBits(len);

			assertThat(readBits, is(i));
		}

		for (int i = Integer.MAX_VALUE - 1000; i < Integer.MAX_VALUE; i++) {
			int len = (int) (1 + Math.log(i) / Math.log(2));
			byte[] buf = Utils.toBytes(i << 64 - len);
			ByteArrayInputStream bais = new ByteArrayInputStream(buf);
			DefaultBitInputStream bis = new DefaultBitInputStream(bais);
			int readBits = bis.readBits(len);

			assertThat(readBits, is(i));
		}
	}

	@Test
	public void test_ReadBits_long_range() throws IOException {
		for (long i = 1L; i < 1000000L; i++) {
			int len = (int) (1 + Math.log(i) / Math.log(2));
			byte[] buf = Utils.toBytes(i << 64 - len);
			ByteArrayInputStream bais = new ByteArrayInputStream(buf);
			DefaultBitInputStream bis = new DefaultBitInputStream(bais);
			long readBits = bis.readLongBits(len);

			assertThat(readBits, is(i));
		}

		for (long i = Integer.MAX_VALUE - 1000L; i < Integer.MAX_VALUE + 1000L; i++) {
			int len = (int) (1 + Math.log(i) / Math.log(2));
			byte[] buf = Utils.toBytes(i << 64 - len);
			ByteArrayInputStream bais = new ByteArrayInputStream(buf);
			DefaultBitInputStream bis = new DefaultBitInputStream(bais);
			long readBits = bis.readLongBits(len);

			assertThat(readBits, is(i));
		}

		for (long i = Long.MAX_VALUE - 1000L; i < Long.MAX_VALUE; i++) {
			int len = (int) (1 + Math.log(i) / Math.log(2));
			byte[] buf = Utils.toBytes(i << 64 - len);
			ByteArrayInputStream bais = new ByteArrayInputStream(buf);
			DefaultBitInputStream bis = new DefaultBitInputStream(bais);
			long readBits = bis.readLongBits(len);

			assertThat(readBits, is(i));
		}
	}
}
