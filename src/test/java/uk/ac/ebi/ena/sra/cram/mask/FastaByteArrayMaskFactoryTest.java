package uk.ac.ebi.ena.sra.cram.mask;

import static org.hamcrest.CoreMatchers.*;
import static org.junit.Assert.assertThat;

import org.junit.Test;

public class FastaByteArrayMaskFactoryTest {
	private FastaByteArrayMaskFactory factory = new FastaByteArrayMaskFactory();

	@Test
	public void testEmptyLine() throws ReadMaskFormatException {
		String line = "";
		PositionMask mask = factory.createMask(line);

		assertThat(mask, notNullValue());
		assertThat(mask.getMaskedPositions(), equalTo(new int[] {}));
	}

	@Test(expected = ReadMaskFormatException.class)
	public void failEmptyLineWithSlashN() throws ReadMaskFormatException {
		String line = "\n";
		factory.createMask(line);
	}

	@Test
	public void test_1() throws ReadMaskFormatException {
		String line = new String(
				new byte[] { FastaByteArrayMaskFactory.DEFAULT_MASK_BYTE });
		PositionMask mask = factory.createMask(line);

		assertThat(mask, notNullValue());
		assertThat(mask.getMaskedPositions(), equalTo(new int[] { 1 }));
	}

	@Test(expected = ReadMaskFormatException.class)
	public void fail_1n() throws ReadMaskFormatException {
		String line = new String(
				new byte[] { FastaByteArrayMaskFactory.DEFAULT_MASK_BYTE, '\n' });
		factory.createMask(line);
	}

	@Test
	public void test_xyxxyyx() throws ReadMaskFormatException {
		String line = "xyxxyyx" ;
		PositionMask mask = factory.createMask(line);

		assertThat(mask, notNullValue());
		assertThat(mask.getMaskedPositions(), equalTo(new int[] { 1, 3, 4, 7 }));
	}

}
