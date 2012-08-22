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
package uk.ac.ebi.ena.sra.cram.impl;

import static org.hamcrest.CoreMatchers.*;
import static org.junit.Assert.assertThat;

import java.io.IOException;
import java.util.ArrayList;

import org.junit.Test;

import uk.ac.ebi.ena.sra.cram.encoding.BaseChange;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;

public class RestoreBasesTest {

	@Test
	public void test1() throws IOException {
		CramRecord record = new CramRecord();
		record.setAlignmentStart(1);
		record.setReadLength(1);
		RestoreBases restorator = new RestoreBases();
		byte[] readBases = "G".getBytes();
		byte[] refBases = "A".getBytes();
		restorator.setProvider(new ByteArraySequenceBaseProvider(refBases));
		SubstitutionVariation v = new SubstitutionVariation();
		v.setPosition(1);
		v.setBaseChange(new BaseChange((byte) 'A', (byte) 'G'));
		record.setReadFeatures(new ArrayList<ReadFeature>());
		record.getReadFeatures().add(v);

		byte[] restoredBases = restorator.restoreReadBases(record);
		assertThat(restoredBases, equalTo(readBases));
	}

	@Test
	public void test2() throws IOException {
		CramRecord record = new CramRecord();
		record.setAlignmentStart(1);
		record.setReadLength(2);
		RestoreBases restorator = new RestoreBases();
		byte[] readBases = "AG".getBytes();
		byte[] refBases = "AA".getBytes();
		restorator.setProvider(new ByteArraySequenceBaseProvider(refBases));
		SubstitutionVariation v = new SubstitutionVariation();
		v.setPosition(2);
		v.setBaseChange(new BaseChange((byte) 'A', (byte) 'G'));
		record.setReadFeatures(new ArrayList<ReadFeature>());
		record.getReadFeatures().add(v);

		byte[] restoredBases = restorator.restoreReadBases(record);
		assertThat(restoredBases, equalTo(readBases));
	}

	@Test
	public void test3() throws IOException {
		CramRecord record = new CramRecord();
		record.setReadFeatures(new ArrayList<ReadFeature>());
		record.setAlignmentStart(1);
		record.setReadLength(5);
		RestoreBases restorator = new RestoreBases();
		byte[] readBases = "AGATA".getBytes();
		byte[] refBases = "AAAAA".getBytes();
		restorator.setProvider(new ByteArraySequenceBaseProvider(refBases));
		SubstitutionVariation v = new SubstitutionVariation();
		v.setPosition(2);
		v.setBaseChange(new BaseChange((byte) 'A', (byte) 'G'));
		record.getReadFeatures().add(v);
		v = new SubstitutionVariation();
		v.setPosition(4);
		v.setBaseChange(new BaseChange((byte) 'A', (byte) 'T'));
		record.getReadFeatures().add(v);

		byte[] restoredBases = restorator.restoreReadBases(record);
		assertThat(restoredBases, equalTo(readBases));
	}

	@Test
	public void test4() throws IOException {
		CramRecord record = new CramRecord();
		record.setReadFeatures(new ArrayList<ReadFeature>());
		record.setAlignmentStart(1);
		record.setReadLength(5);
		RestoreBases restorator = new RestoreBases();
		byte[] readBases = "AGATA".getBytes();
		byte[] refBases = "AAAAA".getBytes();
		restorator.setProvider(new ByteArraySequenceBaseProvider(refBases));
		SubstitutionVariation v = new SubstitutionVariation();
		v.setPosition(2);
		v.setBaseChange(new BaseChange((byte) 'A', (byte) 'G'));
		record.getReadFeatures().add(v);
		v = new SubstitutionVariation();
		v.setPosition(4);
		v.setBaseChange(new BaseChange((byte) 'A', (byte) 'T'));
		record.getReadFeatures().add(v);

		byte[] restoredBases = restorator.restoreReadBases(record);
		assertThat(restoredBases, equalTo(readBases));
	}

	@Test
	public void test5() throws IOException {
		CramRecord record = new CramRecord();
		record.setReadFeatures(new ArrayList<ReadFeature>());
		record.setAlignmentStart(1);
		record.setReadLength(5);
		RestoreBases restorator = new RestoreBases();
		byte[] readBases = "AGATA".getBytes();
		byte[] refBases = "AAAA".getBytes();
		restorator.setProvider(new ByteArraySequenceBaseProvider(refBases));

		SubstitutionVariation v = new SubstitutionVariation();
		v.setPosition(2);
		v.setBaseChange(new BaseChange((byte) 'A', (byte) 'G'));
		record.getReadFeatures().add(v);

		InsertionVariation iv = new InsertionVariation();
		iv.setPosition(4);
		iv.setSequence("T".getBytes());
		record.getReadFeatures().add(iv);

		byte[] restoredBases = restorator.restoreReadBases(record);
		assertThat(restoredBases, equalTo(readBases));
	}

	@Test
	public void test6() throws IOException {
		String read = "ATAGA";
		String reference = "AAAA";
		byte[] readBases = read.getBytes();
		byte[] refBases = reference.getBytes();

		CramRecord record = new CramRecord();
		record.setReadFeatures(new ArrayList<ReadFeature>());
		record.setAlignmentStart(1);
		record.setReadLength(readBases.length);
		RestoreBases restorator = new RestoreBases();
		restorator.setProvider(new ByteArraySequenceBaseProvider(refBases));

		InsertionVariation iv = new InsertionVariation();
		iv.setPosition(2);
		iv.setSequence("T".getBytes());
		record.getReadFeatures().add(iv);

		SubstitutionVariation v = new SubstitutionVariation();
		v.setPosition(4);
		v.setBaseChange(new BaseChange((byte) 'A', (byte) 'G'));
		record.getReadFeatures().add(v);

		byte[] restoredBases = restorator.restoreReadBases(record);
		assertThat(restoredBases, equalTo(readBases));
	}

	@Test
	public void test7() throws IOException {
		String read = "AGT";
		String reference = "ACGT";
		byte[] readBases = read.getBytes();
		byte[] refBases = reference.getBytes();

		CramRecord record = new CramRecord();
		record.setReadFeatures(new ArrayList<ReadFeature>());
		record.setAlignmentStart(1);
		record.setReadLength(readBases.length);

		RestoreBases restorator = new RestoreBases();
		restorator.setProvider(new ByteArraySequenceBaseProvider(refBases));

		DeletionVariation dv = new DeletionVariation();
		dv.setPosition(2);
		dv.setLength(1);
		record.getReadFeatures().add(dv);

		byte[] restoredBases = restorator.restoreReadBases(record);
		assertThat(restoredBases, equalTo(readBases));
	}

	@Test
	public void test8() throws IOException {
		String read = "AGA";
		String reference = "ACGT";
		byte[] readBases = read.getBytes();
		byte[] refBases = reference.getBytes();

		CramRecord record = new CramRecord();
		record.setReadFeatures(new ArrayList<ReadFeature>());
		record.setAlignmentStart(1);
		record.setReadLength(readBases.length);

		RestoreBases restorator = new RestoreBases();
		restorator.setProvider(new ByteArraySequenceBaseProvider(refBases));

		DeletionVariation dv = new DeletionVariation();
		dv.setPosition(2);
		dv.setLength(1);
		record.getReadFeatures().add(dv);

		SubstitutionVariation sv = new SubstitutionVariation();
		sv.setPosition(3);
		sv.setBaseChange(new BaseChange((byte) 'T', (byte) 'A'));
		record.getReadFeatures().add(sv);

		byte[] restoredBases = restorator.restoreReadBases(record);
		assertThat(restoredBases, equalTo(readBases));
	}

	@Test
	public void test9() throws IOException {
		// read: A-GCAGC--CGT
		// ref.: ACGT-CGTA--T

		String readAlignment = "A-GCAGC--CGT";
		String referenceAlignment = "ACGT-CGTA--T";
		assertThat(readAlignment.length(), is(referenceAlignment.length()));

		String read = readAlignment.replaceAll("-", "");
		String reference = referenceAlignment.replaceAll("-", "");
		byte[] readBases = read.getBytes();
		byte[] refBases = reference.getBytes();

		CramRecord record = new CramRecord();
		record.setReadFeatures(new ArrayList<ReadFeature>());
		record.setAlignmentStart(1);
		record.setReadLength(readBases.length);

		RestoreBases restorator = new RestoreBases();
		restorator.setProvider(new ByteArraySequenceBaseProvider(refBases));

		DeletionVariation dv = new DeletionVariation();
		dv.setPosition(2);
		dv.setLength(1);
		record.getReadFeatures().add(dv);

		SubstitutionVariation sv = new SubstitutionVariation();
		sv.setPosition(3);
		sv.setBaseChange(new BaseChange((byte) 'T', (byte) 'C'));
		record.getReadFeatures().add(sv);

		InsertionVariation iv = new InsertionVariation();
		iv.setPosition(4);
		iv.setSequence("A".getBytes());
		record.getReadFeatures().add(iv);

		sv = new SubstitutionVariation();
		sv.setPosition(5);
		sv.setBaseChange(new BaseChange((byte) 'C', (byte) 'G'));
		record.getReadFeatures().add(sv);

		sv = new SubstitutionVariation();
		sv.setPosition(6);
		sv.setBaseChange(new BaseChange((byte) 'G', (byte) 'C'));
		record.getReadFeatures().add(sv);

		dv = new DeletionVariation();
		dv.setPosition(7);
		dv.setLength(2);
		record.getReadFeatures().add(dv);

		iv = new InsertionVariation();
		iv.setPosition(7);
		iv.setSequence("CG".getBytes());
		record.getReadFeatures().add(iv);

		byte[] restoredBases = restorator.restoreReadBases(record);
		assertThat(restoredBases, equalTo(readBases));
	}

	@Test
	public void test10() throws IOException {
		// read: A-GCAGC--CGT
		// ref.: ACGT-CGTA--T

		String readAlignment = "A-GCAGC--CGT";
		String referenceAlignment = "ACGT-CGTA--T";
		assertThat(readAlignment.length(), is(referenceAlignment.length()));

		String read = readAlignment.replaceAll("-", "");
		String reference = referenceAlignment.replaceAll("-", "");
		byte[] readBases = read.getBytes();
		byte[] refBases = reference.getBytes();

		CramRecord record = new CramRecord();
		record.setReadFeatures(new ArrayList<ReadFeature>());
		record.setAlignmentStart(1);
		record.setReadLength(readBases.length);

		RestoreBases restorator = new RestoreBases();
		restorator.setProvider(new ByteArraySequenceBaseProvider(refBases));

		DeletionVariation dv = new DeletionVariation();
		dv.setPosition(2);
		dv.setLength(1);
		record.getReadFeatures().add(dv);

		SubstitutionVariation sv = new SubstitutionVariation();
		sv.setPosition(3);
		sv.setBaseChange(new BaseChange((byte) 'T', (byte) 'C'));
		record.getReadFeatures().add(sv);

		InsertionVariation iv = new InsertionVariation();
		iv.setPosition(4);
		iv.setSequence("A".getBytes());
		record.getReadFeatures().add(iv);

		sv = new SubstitutionVariation();
		sv.setPosition(5);
		sv.setBaseChange(new BaseChange((byte) 'C', (byte) 'G'));
		record.getReadFeatures().add(sv);

		sv = new SubstitutionVariation();
		sv.setPosition(6);
		sv.setBaseChange(new BaseChange((byte) 'G', (byte) 'C'));
		record.getReadFeatures().add(sv);

		dv = new DeletionVariation();
		dv.setPosition(7);
		dv.setLength(2);
		record.getReadFeatures().add(dv);

		iv = new InsertionVariation();
		iv.setPosition(7);
		iv.setSequence("CG".getBytes());
		record.getReadFeatures().add(iv);

		ReverseVariations rv = new ReverseVariations();
		rv.reverse(record);
		rv.reverse(record);

		byte[] restoredBases = restorator.restoreReadBases(record);
		assertThat(restoredBases, equalTo(readBases));
	}

}
