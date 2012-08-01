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

public class RecordPointer extends BitPointer implements
		Comparable<RecordPointer> {
	private long blockStart;
	private long alignmentStart;
	private long recordNumber;

	@Override
	public int compareTo(RecordPointer o) {
		return (int) (alignmentStart - o.alignmentStart);
	}

	public long getBlockStart() {
		return blockStart;
	}

	public void setBlockStart(long blockStart) {
		this.blockStart = blockStart;
	}

	public long getAlignmentStart() {
		return alignmentStart;
	}

	public void setAlignmentStart(long alignmentStart) {
		this.alignmentStart = alignmentStart;
	}

	public long getRecordNumber() {
		return recordNumber;
	}

	public void setRecordNumber(long recordNumber) {
		this.recordNumber = recordNumber;
	}

	@Override
	public String toString() {
		return blockStart + ">" +recordNumber + ":"+ super.toString() + "@" + alignmentStart;
	}
}
