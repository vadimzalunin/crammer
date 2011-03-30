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

	@Test(timeout = 300)
	public void becnhmar_ReadBits_32() throws IOException {
		int maxValues = 1000000;
		byte[] buf = new byte[maxValues * 4];
		byte fillByte = (byte) (64 + 16 + 4 + 1);
		Arrays.fill(buf, (byte) (64 + 16 + 4 + 1));
		int value = 0;
		value = (value << 8) + fillByte;
		value = (value << 8) + fillByte;
		value = (value << 8) + fillByte;
		value = (value << 8) + fillByte;
		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);
		for (int i = 0; i < maxValues; i++)
			bis.readBits(32);

	}

}
