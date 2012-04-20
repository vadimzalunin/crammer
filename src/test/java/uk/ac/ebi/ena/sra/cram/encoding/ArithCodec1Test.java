package uk.ac.ebi.ena.sra.cram.encoding;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.junit.Assert.assertThat;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.Ignore;
import org.junit.Test;

import uk.ac.ebi.ena.sra.cram.format.ByteFrequencies;
import uk.ac.ebi.ena.sra.cram.format.DiByteFrequencies;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.LongBufferBitInputStream;

public class ArithCodec1Test {

	@Test
	@Ignore
	public void test1() throws IOException {
		DiByteFrequencies f = new DiByteFrequencies();
		f.add((byte) 1, (byte) 1, 1);
		f.add((byte) 1, (byte) 2, 2);
		f.add((byte) 2, (byte) 2, 3);

		ArithCodec1 codec = new ArithCodec1(f.getFrequencies(), f.getValues());
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		byte[] bytes = new byte[] { 1, 1, 2, 2 };
		codec.write(bos, bytes);
		bos.flush();
		byte[] buf = baos.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new LongBufferBitInputStream(bais);

		assertThat(codec.read(bis, bytes.length), equalTo(bytes));
	}

	@Test
	@Ignore
	public void test2() throws IOException {
		byte[] message = "qweq;wlekjq;wlkejq;wklj".getBytes();

		DiByteFrequencies f = new DiByteFrequencies();
		for (int i = 1; i < message.length; i++)
			f.add(message[i - 1], message[i]);

		ArithCodec1 codec = new ArithCodec1(f.getFrequencies(), f.getValues());
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		codec.write(bos, message);
		bos.flush();
		byte[] buf = baos.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new LongBufferBitInputStream(bais);

		assertThat(codec.read(bis, message.length), equalTo(message));
	}

	@Test
	public void test3() throws IOException {
		int maxMessages = 1000;
		int minLength = 10;
		int maxLength = 100;
		int alphabetSize = 40;
		int alphabetOffset = 33;

		List<byte[]> messages = new ArrayList<byte[]>(maxMessages);
		DiByteFrequencies f = new DiByteFrequencies();
		ByteFrequencies bf = new ByteFrequencies();

		Random random = new Random();
		for (int m = 0; m < maxMessages; m++) {
			int len = random.nextInt(maxLength - minLength) + minLength;
			byte[] message = new byte[len];
			for (int i = 0; i < len; i++) {
				// message[i] = (byte) (alphabetSize -
				// random.nextInt(alphabetSize) + alphabetOffset);
				message[i] = (byte) (Math.sqrt(alphabetSize * alphabetSize
						- random.nextInt(alphabetSize * alphabetSize)) + alphabetOffset);
				bf.add(message[i]);
			}

			messages.add(message);

			for (int i = 1; i < message.length; i++)
				f.add(message[i - 1], message[i]);
		}
		System.out.println(bf);
		System.out.println(f);

		ArithCodec1 codec = new ArithCodec1(f.getFrequencies(), f.getValues());
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);
		for (int m = 0; m < maxMessages; m++) {
			codec.write(bos, messages.get(m));
		}
		bos.close();
		byte[] buf = baos.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new LongBufferBitInputStream(bais);

		codec = new ArithCodec1(f.getFrequencies(), f.getValues());
		for (int m = 0; m < maxMessages; m++) {
			byte[] message = messages.get(m);
			assertThat(codec.read(bis, message.length), equalTo(message));
			System.out.println(new String(message));
		}
	}

}
