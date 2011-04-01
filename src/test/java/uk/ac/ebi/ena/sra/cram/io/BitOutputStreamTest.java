package uk.ac.ebi.ena.sra.cram.io;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Arrays;

import org.junit.Test;

public class BitOutputStreamTest {

	@Test
	public void test_WriteByte_all_possible_bytes() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);
		for (byte i = Byte.MIN_VALUE; i < Byte.MAX_VALUE; i++) {
			baos.reset();
			bos.writeByte(i);
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
