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

public class CramIndexBuilder {
	private CramIndex index;
	private long alignmentResolution = 10000;

	private int currentSeqIndex = 0;
	private long currentBlockOffset = 0;
	private long recordNumber = 0;

	public CramIndexBuilder(long alignmentResolution) {
		this.alignmentResolution = alignmentResolution;
		index = new CramIndex();
	}

	public void startNewSequence(int seqIndex) {
		currentSeqIndex = seqIndex;
		recordNumber = 0;
	}

	public void startNewBlock(long byteOffset) {
		currentBlockOffset = byteOffset;
	}

	public void addRecord(long alStart, long byteOffset, byte bitOffset) {
		if (recordNumber % alignmentResolution == 0) {
			RecordPointer pointer = new RecordPointer();
			pointer.setAlignmentStart(alStart);
			pointer.setBlockStart(currentBlockOffset);
			pointer.setByteOffset(byteOffset);
			pointer.setBitOffset(bitOffset);
			pointer.setRecordNumber(recordNumber);
			index.addRecordPointer(currentSeqIndex, pointer);
		}
		recordNumber++;
	}

	public CramIndex getCramIndex() {
		return index;
	}

}
