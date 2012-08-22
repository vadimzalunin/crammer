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

import static org.hamcrest.CoreMatchers.*;
import static org.junit.Assert.assertThat;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Arrays;

import org.junit.Test;

import uk.ac.ebi.ena.sra.cram.Utils;

public class BitOutputStreamTest {

	@Test
	public void test_write_long_64bits() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);
		long value = 0;
		value = (1L << 63) | (1L << 62);
		value = value | 1;

		bos.write(value, 64);
		bos.flush();
		byte[] buf = baos.toByteArray();
		assertThat(buf.length, is(8));
		assertThat(buf, equalTo(Utils.toBytes(value)));
		assertThat(Utils.toBitString(buf), equalTo("1100000000000000000000000000000000000000000000000000000000000001"));
	}

	@Test
	public void test_write_long_128_8_bits() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);
		long value = 128;

		bos.write(1, 7);
		bos.write(value, 8);
		bos.flush();
		byte[] buf = baos.toByteArray();
		assertThat(buf.length, is(2));
		assertThat(buf, equalTo(new byte[] { (byte) 3, 0 }));
		assertThat(Utils.toBitString(buf), equalTo("0000001100000000"));
	}

	@Test
	public void test_write_long_10bits() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);
		long value = 0;
		value = (1L << 9) | (1L << 8);
		// System.out.println(Utils.toBitString(Utils.toBytes(value)));
		value = value | 2;
		// System.out.println(Utils.toBitString(Utils.toBytes(value)));

		bos.write(value, 10);
		bos.flush();
		byte[] buf = baos.toByteArray();
		// System.out.println(Utils.toBitString(buf));
		assertThat(buf.length, is(2));
		assertThat(buf, equalTo(new byte[] { -64, -128 }));
		assertThat(Utils.toBitString(buf), equalTo("1100000010000000"));
	}

	@Test
	public void test_write_long_2bits() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);
		long value = 0;
		value = (1L << 63) | (1L << 62);
		value = value | 1;

		bos.write(value, 2);
		bos.flush();
		byte[] buf = baos.toByteArray();
		// System.out.println(Utils.toBitString(buf));
		assertThat(buf.length, is(1));
		assertThat(buf, equalTo(new byte[] { 64 }));
		assertThat(Utils.toBitString(buf), equalTo("01000000"));
	}

	@Test
	public void test_write_long_8bits() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);
		long value = 0;
		value = (1L << 63) | (1L << 62);
		value = value | 1;

		bos.write(value, 8);
		bos.flush();
		byte[] buf = baos.toByteArray();
		// System.out.println(Utils.toBitString(buf));
		assertThat(buf.length, is(1));
		assertThat(buf, equalTo(new byte[] { 1 }));
		assertThat(Utils.toBitString(buf), equalTo("00000001"));
	}

	@Test(timeout = 600)
	public void benchmark_write_long() throws IOException {
		int maxValues = 1000000;
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);
		long value = 0;
		value = (1L << 63) | (1L << 62);
		value = value | 1;

		for (int i = 0; i < maxValues; i++)
			bos.write(value, 64);

		bos.flush();
	}

	@Test
	public void test_write_int_32bits() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);
		int value = 0;
		value = (1 << 31) | (1 << 30);
		value = value | 1;

		bos.write_int_LSB_0(value, 32);
		bos.flush();
		byte[] buf = baos.toByteArray();
		assertThat(buf.length, is(4));
		assertThat(buf, equalTo(Utils.toBytes(value)));
		assertThat(Utils.toBitString(buf), equalTo("11000000000000000000000000000001"));
	}

	@Test
	public void test_write_int_10bits() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);
		int value = 0;
		value = (1 << 9) | (1 << 8);
		// System.out.println(Utils.toBitString(Utils.toBytes(value)));
		value = value | 2;
		// System.out.println(Utils.toBitString(Utils.toBytes(value)));

		bos.write_int_LSB_0(value, 10);
		bos.flush();
		byte[] buf = baos.toByteArray();
		// System.out.println(Utils.toBitString(buf));
		assertThat(buf.length, is(2));
		assertThat(buf, equalTo(new byte[] { -64, -128 }));
		assertThat(Utils.toBitString(buf), equalTo("1100000010000000"));
	}

	@Test
	public void test_write_int_2bits() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);
		int value = 0;
		value = (1 << 31) | (1 << 30);
		value = value | 1;

		bos.write_int_LSB_0(value, 2);
		bos.flush();
		byte[] buf = baos.toByteArray();
		// System.out.println(Utils.toBitString(buf));
		assertThat(buf.length, is(1));
		assertThat(buf, equalTo(new byte[] { 64 }));
		assertThat(Utils.toBitString(buf), equalTo("01000000"));
	}

	@Test
	public void test_write_int_8bits() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);
		int value = 0;
		value = (1 << 31) | (1 << 30);
		value = value | 1;

		bos.write_int_LSB_0(value, 8);
		bos.flush();
		byte[] buf = baos.toByteArray();
		assertThat(buf.length, is(1));
		assertThat(buf, equalTo(new byte[] { 1 }));
		assertThat(Utils.toBitString(buf), equalTo("00000001"));
	}

	@Test(timeout = 300)
	public void benchmark_write_int() throws IOException {
		int maxValues = 1000000;
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);
		int value = 0;
		value = (1 << 31) | (1 << 30);
		value = value | 1;

		for (int i = 0; i < maxValues; i++)
			bos.write_int_LSB_0(value, 32);

		bos.flush();
	}

	@Test
	public void test_write_byte_leftmost_bit() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);
		baos.reset();
		bos.write(1, 1);
		bos.flush();
		assertThat(baos.toByteArray()[0], is((byte) -128));
	}

	@Test
	public void test_WriteByte_all_possible_bytes() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);
		for (byte i = Byte.MIN_VALUE; i < Byte.MAX_VALUE; i++) {
			baos.reset();
			bos.write(i, 8);
			bos.flush();
			assertThat(baos.toByteArray()[0], is(i));
		}
	}

	@Test
	public void test_write_long_higest_bit() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);

		long value = 1 << 64;
		bos.write(value, 1);
		bos.flush();
		byte[] buf = baos.toByteArray();

		assertThat(buf.length, is(1));
		assertThat(buf[0], is((byte) (1 << 7)));
	}

	@Test
	public void test_write_long_lowest_bit() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);

		long value = 1;
		bos.write(value, 64);
		bos.flush();
		byte[] buf = baos.toByteArray();

		assertThat(buf.length, is(8));
		assertThat(buf[7], is((byte) 1));
	}

	@Test
	public void test_WriteBits_2_is_64_due_to_flush() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		byte value = 1 << 1;
		bos.write(value, 3);
		bos.flush();

		assertThat(baos.toByteArray()[0], is((byte) (value << (8 - 3))));
	}

	@Test
	public void test_Write_int_32_bits() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);

		int value = 1 << 1;
		bos.write(value, 32);
		bos.flush();

		byte[] buf = baos.toByteArray();

		assertThat(buf.length, is(4));
		assertThat(buf[0], is((byte) 0));
		assertThat(buf[1], is((byte) 0));
		assertThat(buf[2], is((byte) 0));
		assertThat(buf[3], is((byte) value));
	}

	@Test(timeout = 350)
	public void benchmark_Write_int_32_bits() throws IOException {
		int maxValues = 1000000;
		byte[] buf = new byte[maxValues * 4];
		byte fillByte = (byte) (64 + 16 + 4 + 1);
		Arrays.fill(buf, (byte) (64 + 16 + 4 + 1));
		int value = 0;
		value = (value << 8) + fillByte;
		value = (value << 8) + fillByte;
		value = (value << 8) + fillByte;
		value = (value << 8) + fillByte;

		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);

		for (int i = 0; i < maxValues; i++)
			bos.write(value, 32);

		bos.flush();
	}
}
