package uk.ac.ebi.ena.sra.cram.encoding;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.junit.Assert.assertThat;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import org.junit.Test;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

public class BaseChangeCodecTest {

	@Test
	public void test_BaseChangeConstructorsWorkSimilar() {
		BaseChange change1 = new BaseChange((byte) 'A', (byte) 'C');

		BaseChange change2 = new BaseChange(change1.getChange());

		assertThat(change2, equalTo(change1));
	}

	@Test
	public void test_A_to_C_change() throws IOException {
		BaseChangeCodec codec = new BaseChangeCodec();

		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bitOutputStream = new DefaultBitOutputStream(baos);

		BaseChange change1 = new BaseChange((byte) 'A', (byte) 'C');
		codec.write(bitOutputStream, change1);

		bitOutputStream.flush();

		BitInputStream bitInputStream = new DefaultBitInputStream(
				new ByteArrayInputStream(baos.toByteArray()));
		BaseChange change2 = codec.read(bitInputStream);

		assertThat(change2, equalTo(change1));
	}

	@Test
	public void test_all_possible_substituons() throws IOException {
		BaseChangeCodec codec = new BaseChangeCodec();

		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bitOutputStream = new DefaultBitOutputStream(baos);

		byte[] refBaseArray = "AAAACCCCGGGGTTTTNNNN".getBytes();
		byte[] readBaseArray = "CGTNAGTNACTNACGNACGT".getBytes();

		for (int i = 0; i < refBaseArray.length; i++) {
			BaseChange change = new BaseChange(refBaseArray[i],
					readBaseArray[i]);
			codec.write(bitOutputStream, change);
			System.out.printf("%c\t%c\t%d\n", (char) refBaseArray[i],
					(char) readBaseArray[i], change.getChange());
		}
		bitOutputStream.flush();

		BitInputStream bitInputStream = new DefaultBitInputStream(
				new ByteArrayInputStream(baos.toByteArray()));

		for (int i = 0; i < refBaseArray.length; i++) {
			BaseChange change1 = new BaseChange(refBaseArray[i],
					readBaseArray[i]);
			BaseChange change2 = codec.read(bitInputStream);
			assertThat(change2, equalTo(change1));
		}

	}

}
