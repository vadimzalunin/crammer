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
	public void test1() throws IOException {
		UnaryCodec codec = new UnaryCodec(false, 0);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		long value = 0L;
		long len = codec.write(bos, value);
		assertThat(len, is(value + 1));

		bos.flush();
		byte[] buf = baos.toByteArray();
		String bitString = Utils.toBitString(buf);

		assertThat(bitString, equalTo("10000000"));

		BitInputStream bis = new DefaultBitInputStream(
				new ByteArrayInputStream(buf));
		Long readValue = codec.read(bis);
		assertThat(readValue, equalTo(value));
	}

	@Test
	public void test2() throws IOException {
		UnaryCodec codec = new UnaryCodec(false, 0);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		long value = 20L;
		long len = codec.write(bos, value);
		assertThat(len, is(value + 1));

		bos.flush();
		byte[] buf = baos.toByteArray();
		String bitString = Utils.toBitString(buf);

		assertThat(bitString, equalTo("111111111111111111111000"));

		BitInputStream bis = new DefaultBitInputStream(
				new ByteArrayInputStream(buf));
		Long readValue = codec.read(bis);
		assertThat(readValue, equalTo(value));
	}

	@Test
	public void test3() throws IOException {
		UnaryCodec codec = new UnaryCodec(false, 0);

		for (long i = 0; i < 1000; i++) {
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
			BitOutputStream bos = new DefaultBitOutputStream(baos);
			long len = codec.write(bos, i);
			assertThat(len, is(i + 1));

			bos.flush();
			byte[] buf = baos.toByteArray();
			assertThat((long) buf.length, is(1 + (len >> 3)));

			BitInputStream bis = new DefaultBitInputStream(
					new ByteArrayInputStream(buf));
			Long readValue = codec.read(bis);
			assertThat(readValue, equalTo(i));
		}
	}
	
	@Test
	public void test4() throws IOException {
		long maxValues = 1000000 ;
		
		UnaryCodec codec = new UnaryCodec(false, 0);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		for (long i = 0; i < 1000000; i++) 
			codec.write(bos, 20L);
		
		bos.flush();
		byte[] buf = baos.toByteArray();
		
		BitInputStream bis = new DefaultBitInputStream(
				new ByteArrayInputStream(buf));
		
		
		Long readValue = codec.read(bis);
		assertThat(readValue, equalTo(i));
	}
}
