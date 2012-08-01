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
package uk.ac.ebi.ena.sra.cram.index;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.log4j.Logger;

public class CramIndex {
	private static Logger log = Logger.getLogger(CramIndex.class);

	private Map<String, AlignmentIndex> alignmentIndexes;

	public RecordPointer findRecordPointerAt(String seq, long alStart) {
		AlignmentIndex alignmentIndex = alignmentIndexes.get(seq);
		return alignmentIndex.findRecordPointer(alStart);
	}

	void addRecordPointer(int seqIndex, RecordPointer pointer) {
		AlignmentIndex alignmentIndex = alignmentIndexes.get(seqIndex);
		alignmentIndex.addRecordPointer(pointer);
	}

	public static CramIndex fromFile(File file) throws FileNotFoundException, IOException {
		log.debug("Building cram index from file: " + file.getAbsolutePath());
		InputStream is = new GZIPInputStream(new BufferedInputStream(new FileInputStream(file)));
		return fromTextInputStream(is);
	}

	private static CramIndex fromTextInputStream(InputStream is) {
		long time1 = System.currentTimeMillis();
		Scanner scanner = new Scanner(is);
		CramIndex cramIndex = new CramIndex();
		cramIndex.alignmentIndexes = new TreeMap<String, AlignmentIndex>();
		AlignmentIndex alignmentIndex = null;
		while (scanner.hasNextLine()) {
			String line = scanner.nextLine();
			String[] chunks = line.split("\\s+");

			String seq = chunks[0];
			alignmentIndex = cramIndex.alignmentIndexes.get(seq);
			if (alignmentIndex == null) {
				alignmentIndex = new AlignmentIndex();
				cramIndex.alignmentIndexes.put(seq, alignmentIndex);
			} else
				alignmentIndex = cramIndex.alignmentIndexes.get(seq);

			RecordPointer pointer = new RecordPointer();
			pointer.setBlockStart(Long.valueOf(chunks[2]));
			pointer.setRecordNumber(Long.valueOf(chunks[3]));
			pointer.setByteOffset(Long.valueOf(chunks[4]));
			pointer.setBitOffset(Byte.valueOf(chunks[5]));
			pointer.setAlignmentStart(Long.valueOf(chunks[6]));

			alignmentIndex.addRecordPointer(pointer);
		}

		long time2 = System.currentTimeMillis();
		log.debug("Indexed read in " + (time2 - time1) / 1000 + " seconds.");

		return cramIndex;
	}

	public static void write(CramIndex index, File file) throws FileNotFoundException, IOException {
		OutputStream os = new GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(file)));

		write(index, os);
	}

	private static void write(CramIndex index, OutputStream os) {
		PrintStream ps = new PrintStream(os);
		for (int seqIndex = 0; seqIndex < index.alignmentIndexes.size(); seqIndex++) {
			AlignmentIndex alIndex = index.alignmentIndexes.get(seqIndex);
			for (RecordPointer pointer : alIndex) {
				ps.printf("%d\t%d\t%d\t%d\t%d\t%d\n", seqIndex, pointer.getBlockStart(), pointer.getRecordNumber(),
						pointer.getByteOffset(), pointer.getBitOffset(), pointer.getAlignmentStart());
			}
		}
		ps.flush();
	}
}
