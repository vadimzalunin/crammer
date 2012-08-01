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
package uk.ac.ebi.ena.sra.cram.spot;

import java.io.File;
import java.text.NumberFormat;
import java.util.Iterator;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import org.apache.commons.math.stat.HashMapFrequency;

class TestPairedTemplateAssembler {

	public static void main(String[] args) {
		SAMFileReader
				.setDefaultValidationStringency(ValidationStringency.SILENT);
		SAMFileReader reader = new SAMFileReader(new File(
				"y:/Data/SangerExample/paired/5120_1.bam"));

		PairedTemplateAssembler ta = new PairedTemplateAssembler(1000000,
				100000);

		int maxRecords = 1000000000;
		SAMRecordIterator iterator = reader.query(reader.getFileHeader().getSequence(0).getSequenceName(), -1, -1, true) ;
		HashMapFrequency frequency = new HashMapFrequency();
		HashMapFrequency insertSizeFrequency = new HashMapFrequency();
		while (iterator.hasNext()) {
			SAMRecord samRecord = iterator.next();
			// if
			// (!"IL24_5120:1:29:15064:14308#0".equals(samRecord.getReadName()))
			// continue;

			if (!samRecord.getProperPairFlag()
					|| samRecord.getReadUnmappedFlag()
					|| samRecord.getMateUnmappedFlag())
				continue;

			int insertSize = samRecord.getInferredInsertSize();
			if (insertSize > -1)
				insertSizeFrequency.addValue(insertSize);
			
			ta.addSAMRecord(samRecord);

			while ((samRecord = ta.nextSAMRecord()) != null) {
				int distance = ta.distanceToNextFragment();
				if (distance > -2)
					frequency.addValue(distance);
				// System.out.printf("%s\t%d\t->\t%d\t%d\n",
				// samRecord.getReadName(), samRecord.getAlignmentStart(),
				// samRecord.getMateAlignmentStart(),
				// distance);
			}

			if (maxRecords-- <= 0)
				break;
		}

		// SAMRecord samRecord = null;
		// while ((samRecord = ta.fetchNextSAMRecord()) != null) {
		// samRecord = ta.nextSAMRecord();
		// int distance = ta.distanceToNextFragment();
		// frequency.addValue(distance);
		// // System.out.printf("%s\t%d->\t%d\n", samRecord.getReadName(),
		// // samRecord.getAlignmentStart(),
		// // samRecord.getMateAlignmentStart());
		// }

		System.out.println(toString(frequency, 100000, 1000));
		System.out.println(toString(insertSizeFrequency, 100000, 1));
	}

	private static String toString(HashMapFrequency freqTable, int cutOff, int step) {
		NumberFormat nf = NumberFormat.getPercentInstance();
		StringBuilder outBuffer = new StringBuilder();
		outBuffer.append("Value \t Freq. \t Pct. \t Cum Pct. \n");
		Iterator<Comparable<?>> iter = freqTable.valuesIterator();
		int counter = 0;
		while (iter.hasNext()) {
			Comparable<?> value = iter.next();
			if (counter++ % step != 0)
				continue;
			outBuffer.append(value);
			outBuffer.append('\t');
			outBuffer.append(freqTable.getCount(value));
			outBuffer.append('\t');
			outBuffer.append(nf.format(freqTable.getPct(value)));
			outBuffer.append('\t');
			outBuffer.append(nf.format(freqTable.getCumPct(value)));
			outBuffer.append('\n');
			if (counter++ > cutOff)
				break;
		}
		return outBuffer.toString();
	}
}
