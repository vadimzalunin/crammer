package uk.ac.ebi.ena.sra.cram.impl;

import static org.junit.Assert.assertTrue;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.Test;

import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

public class CramRecordWriteReadTests {

	@Test
	public void test1() throws IOException {
		ByteArrayOutputStream os = new ByteArrayOutputStream(1024);
		BitOutputStream bos = new DefaultBitOutputStream(os);
		CramRecordWriter writer = new CramRecordWriter(bos);

		CramRecord record = new CramRecord();
		record.setAlignmentStart(50);
		record.setNegativeStrand(true);
		record.setPerfectMatch(true);

		writer.writeCramRecord(record);
		writer.flush();

		byte[] data = os.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(data);
		CramRecordIterator reader = new CramRecordIterator(
				new DefaultBitInputStream(bais),
				new ByteArraySequenceBaseProvider(
						"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".getBytes()));
		assertTrue(reader.hasNext());

		CramRecord record2 = reader.next();

		assertTrue(record.toString() + "\n" + record2.toString(),
				record.equals(record2));
	}

	@Test
	public void test2() throws IOException {
		ByteArrayOutputStream os = new ByteArrayOutputStream(1024);
		BitOutputStream bos = new DefaultBitOutputStream(os);
		CramRecordWriter writer = new CramRecordWriter(bos);

		CramRecord record1 = new CramRecord();
		record1.setAlignmentStart(13);
		record1.setNegativeStrand(false);
		record1.setPerfectMatch(true);

		writer.writeCramRecord(record1);

		CramRecord record2 = new CramRecord();
		record2.setAlignmentStart(21);
		record2.setNegativeStrand(false);
		record2.setPerfectMatch(true);

		writer.writeCramRecord(record2);
		writer.flush();

		byte[] data = os.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(data);
		CramRecordIterator reader = new CramRecordIterator(
				new DefaultBitInputStream(bais),
				new ByteArraySequenceBaseProvider(
						"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".getBytes()));
		assertTrue(reader.hasNext());
		CramRecord record3 = reader.next();
		assertTrue(record3.toString() + "\n" + record1.toString(),
				record3.equals(record1));

		assertTrue(reader.hasNext());
		CramRecord record4 = reader.next();
		assertTrue(record4.toString() + "\n" + record2.toString(),
				record4.equals(record2));
	}

	@Test
	public void test3() throws IOException {
		ByteArrayOutputStream os = new ByteArrayOutputStream(1024);
		BitOutputStream bos = new DefaultBitOutputStream(os);
		CramRecordWriter writer = new CramRecordWriter(bos);

		CramRecord record = new CramRecord();
		record.setAlignmentStart(1);
		record.setNegativeStrand(false);
		record.setPerfectMatch(false);

		DeletionVariation del = new DeletionVariation();
		del.setPosition(2);
		del.setLength(1);
		List<DeletionVariation> dels = new ArrayList<DeletionVariation>();
		dels.add(del);
		record.setDeletionVariations(dels);

		writer.writeCramRecord(record);
		writer.flush();

		byte[] data = os.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(data);
		CramRecordIterator reader = new CramRecordIterator(
				new DefaultBitInputStream(bais),
				new ByteArraySequenceBaseProvider(
						"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".getBytes()));
		assertTrue(reader.hasNext());

		CramRecord record2 = reader.next();

		assertTrue(record.toString() + "\n" + record2.toString(),
				record.equals(record2));
	}

	@Test
	public void test4() throws IOException {
		ByteArrayOutputStream os = new ByteArrayOutputStream(1024);
		BitOutputStream bos = new DefaultBitOutputStream(os);
		CramRecordWriter writer = new CramRecordWriter(bos);

		CramRecord record = new CramRecord();
		record.setAlignmentStart(1);
		record.setNegativeStrand(false);
		record.setPerfectMatch(false);

		SubstitutionVariation sub = new SubstitutionVariation();
		sub.setPosition(1);
		sub.setBase((byte) 'C');
		sub.setRefernceBase((byte) 'A');
		List<SubstitutionVariation> subs = new ArrayList<SubstitutionVariation>();
		subs.add(sub);
		record.setSubstitutionVariations(subs);

		writer.writeCramRecord(record);
		writer.flush();

		byte[] data = os.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(data);
		CramRecordIterator reader = new CramRecordIterator(
				new DefaultBitInputStream(bais),
				new ByteArraySequenceBaseProvider(
						"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".getBytes()));
		assertTrue(reader.hasNext());

		CramRecord record2 = reader.next();

		assertTrue(record.toString() + "\n" + record2.toString(),
				record.equals(record2));
	}

	@Test
	public void test5() throws IOException {
		ByteArrayOutputStream os = new ByteArrayOutputStream(1024);
		BitOutputStream bos = new DefaultBitOutputStream(os);
		CramRecordWriter writer = new CramRecordWriter(bos);

		CramRecord record = new CramRecord();
		record.setAlignmentStart(1);
		record.setNegativeStrand(false);
		record.setPerfectMatch(false);

		InsertionVariation ins = new InsertionVariation();
		ins.setPosition(29);
		ins.setSequence("C".getBytes());
		List<InsertionVariation> subs = new ArrayList<InsertionVariation>();
		subs.add(ins);
		record.setInsertionVariations(subs);

		writer.writeCramRecord(record);
		writer.flush();

		byte[] data = os.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(data);
		CramRecordIterator reader = new CramRecordIterator(
				new DefaultBitInputStream(bais),
				new ByteArraySequenceBaseProvider(
						"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".getBytes()));
		assertTrue(reader.hasNext());

		CramRecord record2 = reader.next();

		assertTrue(record.toString() + "\n" + record2.toString(),
				record.equals(record2));
	}

	@Test
	public void test6() throws IOException {
		int maxTests = 1000;
		ByteArrayOutputStream os = new ByteArrayOutputStream(1024);
		BitOutputStream bos = new DefaultBitOutputStream(os);
		CramRecordWriter writer = new CramRecordWriter(bos);

		SequenceBaseProvider referenceBaseProvider = new SequenceBaseProvider() {

			@Override
			public byte getBaseAt(String seqName, long position) {
				return "ACGT".getBytes()[(int) (position % 4)];
			}
		};

		List<CramRecord> records = new ArrayList<CramRecord>();
		Random random = new Random();
		CramRecord prevRecord = null;
		for (int i = 0; i < maxTests; i++) {
			CramRecord record = new CramRecord();
			int randomOffset = 1 + random.nextInt(256);
			if (prevRecord == null)
				record.setAlignmentStart(randomOffset);
			else
				record.setAlignmentStart(prevRecord.getAlignmentStart()
						+ randomOffset);
			prevRecord = record;

			record.setNegativeStrand(random.nextBoolean());
			if (!random.nextBoolean()) {
				record.setPerfectMatch(false);

				do {
					byte[] alphabet = "ACGTN".getBytes();
					byte[] sequence = new byte[1];
					for (int j = 0; j < sequence.length; j++)
						sequence[j] = alphabet[random.nextInt(alphabet.length)];

					InsertionVariation ins = new InsertionVariation();
					ins.setPosition(1 + random.nextInt(31));
					ins.setSequence(sequence);
					List<InsertionVariation> subs = new ArrayList<InsertionVariation>();
					subs.add(ins);
					record.setInsertionVariations(subs);
				} while (random.nextFloat() > 0.9f);
			} else
				record.setPerfectMatch(true);

			writer.writeCramRecord(record);
			records.add(record);
		}

		writer.flush();

		byte[] data = os.toByteArray();
		ByteArrayInputStream bais = new ByteArrayInputStream(data);

		CramRecordIterator reader = new CramRecordIterator(
				new DefaultBitInputStream(bais), referenceBaseProvider);

		for (int i = 0; i < maxTests; i++) {
			assertTrue(reader.hasNext());
			CramRecord record = records.get(i);
			CramRecord record2 = reader.next();

			assertTrue(
					"Test: " + record.toString() + "\n" + record2.toString(),
					records.get(i).equals(record2));
		}

	}

	@Test
	public void test7() throws IOException {
		ByteArrayOutputStream os = new ByteArrayOutputStream(1024);
		BitOutputStream bos = new DefaultBitOutputStream(os);
		CramRecordWriter writer = new CramRecordWriter(bos);

		CramRecord record1 = new CramRecord();
		record1.setAlignmentStart(21);
		record1.setNegativeStrand(false);
		record1.setPerfectMatch(false);
		InsertionVariation ins = new InsertionVariation();
		ins.setPosition(29);
		ins.setSequence("CC".getBytes());
		List<InsertionVariation> subs = new ArrayList<InsertionVariation>();
		subs.add(ins);
		record1.setInsertionVariations(subs);
		writer.writeCramRecord(record1);

		CramRecord record2 = new CramRecord();
		record2.setAlignmentStart(123);
		record2.setNegativeStrand(false);
		record2.setPerfectMatch(true);
		writer.writeCramRecord(record2);

		writer.flush();

		byte[] data = os.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(data);
		CramRecordIterator reader = new CramRecordIterator(
				new DefaultBitInputStream(bais),
				new ByteArraySequenceBaseProvider(
						"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".getBytes()));
		assertTrue(reader.hasNext());
		CramRecord record3 = reader.next();
		assertTrue(record3.toString() + "\n" + record1.toString(),
				record3.equals(record1));

		assertTrue(reader.hasNext());
		CramRecord record4 = reader.next();
		assertTrue(record4.toString() + "\n" + record2.toString(),
				record4.equals(record2));
	}

	@Test
	public void testMultiVars() throws IOException {
		ByteArrayOutputStream os = new ByteArrayOutputStream(1024);
		BitOutputStream bos = new DefaultBitOutputStream(os);
		CramRecordWriter writer = new CramRecordWriter(bos);

		CramRecord record = new CramRecord();
		record.setAlignmentStart(1);
		record.setNegativeStrand(false);
		record.setPerfectMatch(false);

		List<SubstitutionVariation> subs = new ArrayList<SubstitutionVariation>();
		int inReadPos = 1;

		SubstitutionVariation sub = new SubstitutionVariation();
		sub.setPosition(inReadPos);
		sub.setBase((byte) 'C');
		sub.setRefernceBase((byte) 'A');
		subs.add(sub);

		inReadPos += 2;
		sub = new SubstitutionVariation();
		sub.setPosition(inReadPos);
		sub.setBase((byte) 'G');
		sub.setRefernceBase((byte) 'A');
		subs.add(sub);

		record.setSubstitutionVariations(subs);

		List<InsertionVariation> ins = new ArrayList<InsertionVariation>();

		inReadPos += 3;
		InsertionVariation iv = new InsertionVariation();
		iv.setPosition(inReadPos);
		iv.setSequence("CC".getBytes());
		ins.add(iv);

		record.setInsertionVariations(ins);

		writer.writeCramRecord(record);

		CramRecord record2 = new CramRecord();
		record2.setAlignmentStart(record.getAlignmentStart() + 12);
		record2.setNegativeStrand(false);
		record2.setPerfectMatch(true);
		writer.writeCramRecord(record2);

		writer.flush();

		byte[] data = os.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(data);
		CramRecordIterator reader = new CramRecordIterator(
				new DefaultBitInputStream(bais),
				new ByteArraySequenceBaseProvider(
						"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".getBytes()));
		assertTrue(reader.hasNext());

		CramRecord record3 = reader.next();

		assertTrue(record.toString() + "\n" + record3.toString(),
				record.equals(record3));

		assertTrue(reader.hasNext());
		CramRecord record4 = reader.next();

		assertTrue(record2.toString() + "\n" + record4.toString(),
				record2.equals(record4));
	}

	@Test
	public void testSubAndIns() throws IOException {
		ByteArrayOutputStream os = new ByteArrayOutputStream(1024);
		BitOutputStream bos = new DefaultBitOutputStream(os);
		CramRecordWriter writer = new CramRecordWriter(bos);

		CramRecord record1 = new CramRecord();
		record1.setAlignmentStart(1);
		record1.setNegativeStrand(false);
		record1.setPerfectMatch(false);
		InsertionVariation iv = new InsertionVariation();
		iv.setPosition(30);
		iv.setSequence("G".getBytes());
		List<InsertionVariation> ins = new ArrayList<InsertionVariation>();
		ins.add(iv);
		record1.setInsertionVariations(ins);

		SubstitutionVariation sv = new SubstitutionVariation();
		sv.setPosition(34);
		sv.setBase((byte) 'G');
		sv.setRefernceBase((byte) 'A');
		List<SubstitutionVariation> subs = new ArrayList<SubstitutionVariation>();
		subs.add(sv);
		record1.setSubstitutionVariations(subs);

		writer.writeCramRecord(record1);

		CramRecord record2 = new CramRecord();
		record2.setAlignmentStart(2);
		record2.setNegativeStrand(false);
		record2.setPerfectMatch(true);
		writer.writeCramRecord(record2);

		writer.flush();

		byte[] data = os.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(data);
		CramRecordIterator reader = new CramRecordIterator(
				new DefaultBitInputStream(bais),
				new ByteArraySequenceBaseProvider(
						"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACC".getBytes()));
		assertTrue(reader.hasNext());
		CramRecord record3 = reader.next();
		assertTrue(record3.toString() + "\n" + record1.toString(),
				record3.equals(record1));

		assertTrue(reader.hasNext());
		CramRecord record4 = reader.next();
		assertTrue(record4.toString() + "\n" + record2.toString(),
				record4.equals(record2));
	}

	@Test(timeout = 1500)
	public void becnhmark_Write() throws IOException {
		int maxRecords = 1000000;

		ByteArrayOutputStream os = new ByteArrayOutputStream(2 * 1024 * 1024);
		BitOutputStream bos = new DefaultBitOutputStream(os);
		CramRecordWriter writer = new CramRecordWriter(bos);

		CramRecord record1 = new CramRecord();
		record1.setAlignmentStart(1);
		record1.setNegativeStrand(false);
		record1.setPerfectMatch(false);

		SubstitutionVariation sv = new SubstitutionVariation();
		sv.setPosition(10);
		sv.setBase((byte) 'G');
		sv.setRefernceBase((byte) 'A');
		List<SubstitutionVariation> subs = new ArrayList<SubstitutionVariation>();
		subs.add(sv);
		record1.setSubstitutionVariations(subs);
		for (int i = 0; i < maxRecords; i++)
			writer.writeCramRecord(record1);

		writer.flush();
	}

	@Test(timeout = 2000)
	public void benchmark_Read() throws IOException {
		int maxRecords = 1000000;

		ByteArrayOutputStream os = new ByteArrayOutputStream(2 * 1024 * 1024);
		BitOutputStream bos = new DefaultBitOutputStream(os);
		CramRecordWriter writer = new CramRecordWriter(bos);

		CramRecord record1 = new CramRecord();
		record1.setAlignmentStart(1);
		record1.setNegativeStrand(false);
		record1.setPerfectMatch(false);

		SubstitutionVariation sv = new SubstitutionVariation();
		sv.setPosition(10);
		sv.setBase((byte) 'G');
		sv.setRefernceBase((byte) 'A');
		List<SubstitutionVariation> subs = new ArrayList<SubstitutionVariation>();
		subs.add(sv);
		record1.setSubstitutionVariations(subs);
		for (int i = 0; i < maxRecords; i++)
			writer.writeCramRecord(record1);

		writer.flush();

		byte[] data = os.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(data);
		CramRecordIterator reader = new CramRecordIterator(
				new DefaultBitInputStream(bais),
				new ByteArraySequenceBaseProvider(
						"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACC".getBytes()));
		for (int i = 0; i < maxRecords; i++) {
			assertTrue(reader.hasNext());
			reader.next();
		}
	}
}
