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

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.util.Arrays;

import net.sf.picard.reference.ReferenceSequenceFileFactory;
import uk.ac.ebi.ena.sra.cram.CramIndexer;
import uk.ac.ebi.ena.sra.cram.CramTools;

public class TestSAMFileReader {

	public static void main(String[] args) throws Exception {
		File bamFile = new File("data/set5/input.bam");
		File refFile = new File("data/set5/ref.fa");

		File cramFile = File.createTempFile(bamFile.getName(), ".cram");
		cramFile.deleteOnExit();
		File indexFile = File.createTempFile(cramFile.getName(), ".crai");
		indexFile.deleteOnExit();

		String command = String.format(
				"-l error cram --input-bam-file %s --reference-fasta-file %s --output-cram-file %s",
				bamFile.getAbsolutePath(), refFile.getAbsolutePath(), cramFile.getAbsolutePath());
		CramTools.main(command.split("\\s+"));

		indexFile = File.createTempFile(cramFile.getName(), ".crai");
		indexFile.deleteOnExit();
		command = String.format("--input-cram-file %s --reference-fasta-file %s --index-file %s",
				cramFile.getAbsolutePath(), refFile.getAbsolutePath(), indexFile.getAbsolutePath());
		CramIndexer.main(command.split("\\s+"));

		DataInputStream fis = new DataInputStream(new FileInputStream(cramFile));
		byte[] buf = new byte[10];
		fis.readFully(buf);
		System.out.println(new String(buf));
		System.out.println(Arrays.toString(buf));
		fis.close();

		ReferenceDiscovery.referenceFactory.put(cramFile,
				ReferenceSequenceFileFactory.getReferenceSequenceFile(refFile));
		SAMFileReader reader = new SAMFileReader(cramFile, indexFile);
		SAMRecordIterator iterator = reader.iterator();

		long counter = 0;
		while (iterator.hasNext()) {
			SAMRecord record = iterator.next();
			if (counter < 10)
				System.out.println(record.format());
			counter++;
		}

		System.out.println(counter);

		reader = new SAMFileReader(cramFile, indexFile);
		iterator = reader.query("20", 60418, 60473, false);

		counter = 0;
		while (iterator.hasNext()) {
			SAMRecord record = iterator.next();
			System.out.println(record.format());
			counter++;
		}

		System.out.println(counter);
	}
}
