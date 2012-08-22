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

import java.util.Comparator;
import java.util.PriorityQueue;

import net.sf.samtools.SAMRecord;

public class LocalReorderingSAMRecordQueue {
	private PriorityQueue<SAMRecord> delegateQueue;
	private boolean drain = false;
	private int minAlignmentDelay = 100;
	private int maxAlignmentStart = 0;
	private int counter = 0;

	public LocalReorderingSAMRecordQueue(int minAlignmentDelay) {
		this.minAlignmentDelay = minAlignmentDelay;
		delegateQueue = new PriorityQueue<SAMRecord>(1000, alignmentStartComparator);
	}

	public void add(SAMRecord record) {
		record.setAttribute("X~", counter++);
		delegateQueue.add(record);
		maxAlignmentStart = Math.max(maxAlignmentStart, record.getAlignmentStart());
	}

	private SAMRecord doPoll() {
		SAMRecord record = delegateQueue.poll();
		if (record == null)
			return null;
		record.setAttribute("X~", null);
		return record;
	}

	public SAMRecord poll() {
		if (drain)
			return doPoll();

		SAMRecord head = delegateQueue.peek();
		if (head == null)
			return null;

		if (maxAlignmentStart - head.getAlignmentStart() < minAlignmentDelay)
			return null;

		return doPoll();
	}

	public void dump() {
		while (!delegateQueue.isEmpty()) {
			SAMRecord r = doPoll();
			System.out.println(r.getAlignmentStart() + "\t" + r.getReadName());
		}
	}

	public boolean isDrain() {
		return drain;
	}

	public void setDrain(boolean drain) {
		this.drain = drain;
	}

	private static Comparator<SAMRecord> alignmentStartComparator = new Comparator<SAMRecord>() {

		/**
		 * The sorting by alignment start is not stable. Therefore we use an
		 * integere value here, which is incremented for each record added/
		 * 
		 * @param o1
		 * @param o2
		 * @return
		 */
		@Override
		public int compare(SAMRecord o1, SAMRecord o2) {

			int result = o1.getAlignmentStart() - o2.getAlignmentStart();
			if (result != 0)
				return result;

			int c1 = (Integer) o1.getAttribute("X~");
			int c2 = (Integer) o2.getAttribute("X~");
			return c1 - c2;

			// result = o1.getReadLength() - o2.getReadLength();
			// if (result != 0)
			// return result;
			//
			// result = o1.getReadString().compareTo(o2.getReadString());
			// if (result != 0)
			// return result;
			//
			// result = o1.getReadName().compareTo(o2.getReadName());
			// return result;

		}
	};
}
