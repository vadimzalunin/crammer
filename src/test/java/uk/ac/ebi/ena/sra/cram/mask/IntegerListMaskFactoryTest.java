package uk.ac.ebi.ena.sra.cram.mask;

import static org.hamcrest.CoreMatchers.*;
import static org.junit.Assert.assertThat;

import org.junit.Test;

public class IntegerListMaskFactoryTest {

	@Test
	public void testEmptyLine() throws ReadMaskFormatException {
		String line = "";
		IntegerListMaskFactory factory = new IntegerListMaskFactory();
		PositionMask mask = factory.createMask(line);

		assertThat(mask, notNullValue());
		assertThat(mask.getMaskedPositions(), equalTo(new int[] {}));
	}

	@Test(expected = ReadMaskFormatException.class)
	public void failEmptyLineWithSlashN() throws ReadMaskFormatException {
		String line = "\n";
		IntegerListMaskFactory factory = new IntegerListMaskFactory();
		factory.createMask(line);
	}
	
	@Test
	public void test_1() throws ReadMaskFormatException {
		String line = "1";
		IntegerListMaskFactory factory = new IntegerListMaskFactory();
		PositionMask mask = factory.createMask(line);

		assertThat(mask, notNullValue());
		assertThat(mask.getMaskedPositions(), equalTo(new int[] {1}));
	}
	
	@Test(expected = ReadMaskFormatException.class)
	public void fail_1n() throws ReadMaskFormatException {
		String line = "1\n";
		IntegerListMaskFactory factory = new IntegerListMaskFactory();
		factory.createMask(line);
	}
	
	@Test(expected = ReadMaskFormatException.class)
	public void fail_non_digit() throws ReadMaskFormatException {
		String line = "A";
		IntegerListMaskFactory factory = new IntegerListMaskFactory();
		factory.createMask(line);
	}
	
	@Test
	public void test_1_2() throws ReadMaskFormatException {
		String line = "1 2";
		IntegerListMaskFactory factory = new IntegerListMaskFactory();
		PositionMask mask = factory.createMask(line);

		assertThat(mask, notNullValue());
		assertThat(mask.getMaskedPositions(), equalTo(new int[] {1, 2}));
	}
	
	@Test
	public void test_1_3_15() throws ReadMaskFormatException {
		String line = "1 3 15";
		IntegerListMaskFactory factory = new IntegerListMaskFactory();
		PositionMask mask = factory.createMask(line);

		assertThat(mask, notNullValue());
		assertThat(mask.getMaskedPositions(), equalTo(new int[] {1, 3, 15}));
	}

}
