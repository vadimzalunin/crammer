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

public class SubexpCodecTest {

	@Test
	public void test_numberOfBits() {
		SubexpCodec codec = new SubexpCodec(2);
		assertThat(codec.numberOfBits(10L), is(6L));
		codec = new SubexpCodec(3);
		assertThat(codec.numberOfBits(14L), is(5L));
	}
	
	@Test
	public void test_write_0_3() throws IOException {
		SubexpCodec codec = new SubexpCodec(3);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);
		long value = 0L;
		int len = (int) codec.write(bos, value);
		assertThat(len, is(4));
		bos.flush();

		byte[] buf = baos.toByteArray();
		assertThat(Utils.toBitString(buf).substring(0, len), equalTo("0000"));
		assertThat(Utils.toBitString(buf), equalTo("00000000"));
	}
	
	@Test
	public void test_write_1_3() throws IOException {
		SubexpCodec codec = new SubexpCodec(3);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);
		long value = 1L;
		int len = (int) codec.write(bos, value);
		assertThat(len, is(4));
		bos.flush();

		byte[] buf = baos.toByteArray();
		assertThat(Utils.toBitString(buf).substring(0, len), equalTo("0001"));
		assertThat(Utils.toBitString(buf), equalTo("00010000"));
	}
	
	@Test
	public void test_write_14_3() throws IOException {
		SubexpCodec codec = new SubexpCodec(3);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);
		long value = 14L;
		int len = (int) codec.write(bos, value);
		assertThat(len, is(5));
		bos.flush();

		byte[] buf = baos.toByteArray();
		assertThat(Utils.toBitString(buf).substring(0, len), equalTo("10110"));
		assertThat(Utils.toBitString(buf), equalTo("10110000"));
	}
	
	@Test
	public void test_read_14_3() throws IOException {
		int m = 3;
		long value = 14L;
		// 10110:
		byte[] buf = new byte[]{(byte) (22<<3)} ;
		
		SubexpCodec codec = new SubexpCodec(m);
		ByteArrayInputStream bais = new ByteArrayInputStream(buf) ;
		BitInputStream bis = new DefaultBitInputStream(bais) ;
		
		
		Long readValue = codec.read(bis);
		assertThat(readValue, is(value));
	}
	
	@Test
	public void test_write_14_4() throws IOException {
		SubexpCodec codec = new SubexpCodec(4);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);
		long value = 14L;
		int len = (int) codec.write(bos, value);
		assertThat(len, is(5));
		bos.flush();

		byte[] buf = baos.toByteArray();
		assertThat(Utils.toBitString(buf).substring(0, len), equalTo("01110"));
		assertThat(Utils.toBitString(buf), equalTo("01110000"));
	}
	
	@Test
	public void test_read_14_4() throws IOException {
		int m = 4;
		long value = 14L;
		// 01110:
		byte[] buf = new byte[]{(byte) (14<<3)} ;
		
		SubexpCodec codec = new SubexpCodec(m);
		ByteArrayInputStream bais = new ByteArrayInputStream(buf) ;
		BitInputStream bis = new DefaultBitInputStream(bais) ;
		
		
		Long readValue = codec.read(bis);
		assertThat(readValue, is(value));
	}

}
