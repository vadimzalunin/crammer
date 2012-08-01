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
package net.sf.samtools;

import java.io.File;
import java.io.FileInputStream;
import java.util.Arrays;

import uk.ac.ebi.ena.sra.cram.CramIndexer;
import uk.ac.ebi.ena.sra.cram.Utils;

public class TestCramPicardIntegration {

	public static void main(String[] args) throws Exception {
		File inputBamFile = new File(args[0]);

		SAMFileReader reader = new SAMFileReader(inputBamFile);
		SAMFileHeader header = reader.getFileHeader();

		File cramOutputFile = File.createTempFile("output", ".cram");
		SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, cramOutputFile);

		copy(reader, writer);

		writer.close();
		reader.close();

		FileInputStream fis = new FileInputStream(cramOutputFile);
		if (!Utils.isCRAM(fis))
			throw new RuntimeException("Written file is not in CRAM format.");
		fis.close();

		CramIndexer.index(cramOutputFile, 1000);

		SAMFileReader bamReader = new SAMFileReader(inputBamFile);
		SAMFileReader cramReader = new SAMFileReader(cramOutputFile);

		SAMRecordIterator bamIterator = bamReader.iterator();
		SAMRecordIterator cramIterator = cramReader.iterator();

		testAlignmentEquality(bamIterator, cramIterator);

		bamReader.close();
		cramReader.close();

		testIndex(inputBamFile, cramOutputFile);

	}

	private static void copy(SAMFileReader reader, SAMFileWriter writer) {
		SAMRecordIterator iterator = reader.iterator();
		while (iterator.hasNext())
			writer.addAlignment(iterator.next());
	}

	private static void testAlignmentEquality(SAMRecordIterator iterator1, SAMRecordIterator iterator2) {
		long counter = 0;
		while (iterator1.hasNext()) {
			SAMRecord bamRecord = iterator1.next();

			if (!iterator2.hasNext()) {
				System.err.println(bamRecord.getSAMString());
				throw new RuntimeException("Record not found in CRAM.");
			}

			SAMRecord cramRecord = iterator2.next();
			
			if (bamRecord.getAlignmentStart() != cramRecord.getAlignmentStart()) {
				System.err.println("Broken at record "+counter);
				System.err.println(bamRecord.getSAMString());
				System.err.println(cramRecord.getSAMString());
				throw new RuntimeException("Alignment start does not match.");
			}
			
			if (!Arrays.equals(bamRecord.getReadBases(), cramRecord.getReadBases())) {
				System.err.println("Broken at record "+counter);
				System.err.println(bamRecord.getSAMString());
				System.err.println(cramRecord.getSAMString());
				throw new RuntimeException("Bases don't match.");
			}
			
			counter++;
//			System.err.printf("%d\t%s\n", bamRecord.getAlignmentStart(), new String (bamRecord.getReadBases()));
		}
		if (iterator2.hasNext())
			throw new RuntimeException("Extra records in iterator 2.");

		System.err.println("Tested " + counter + " records ok.");
	}

	private static void testIndex(File bamFile, File cramFile) {
		SAMFileReader bamReader = new SAMFileReader(bamFile, new File(bamFile.getAbsolutePath() + ".bai"));
		SAMRecordIterator iterator = bamReader.iterator();
		SAMRecord firstSamRecord = iterator.next();
		String referenceName = firstSamRecord.getReferenceName();
		int counter = 100000;
		long firstAlStart = firstSamRecord.getAlignmentStart();
		long lastAlStart = -1;
		while (iterator.hasNext() && counter-- > 0) {
			SAMRecord samRecord = iterator.next();
			if (!referenceName.equals(samRecord.getReferenceName()))
				break;
			lastAlStart = samRecord.getAlignmentStart();
		}
		iterator.close();

		int from = (int) ((lastAlStart - firstAlStart) / 2);
		int to = (int) ((lastAlStart - firstAlStart) / 3 * 2);
		System.err.println("Query range: " + from + " - " + to);

		SAMRecordIterator bamQuery = bamReader.query(referenceName, from, to, false);

		SAMFileReader cramReader = new SAMFileReader(cramFile, new File(cramFile.getAbsolutePath() + ".crai"));
		SAMRecordIterator cramQuery = cramReader.query(referenceName, from, to, false);

		testAlignmentEquality(bamQuery, cramQuery);

		bamQuery.close();
		cramQuery.close();

		bamReader.close();
		cramReader.close();
	}
}
