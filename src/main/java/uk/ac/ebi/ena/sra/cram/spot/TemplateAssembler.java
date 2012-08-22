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

import java.util.Comparator;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.TreeMap;

import net.sf.samtools.SAMRecord;

class TemplateAssembler {
	private int alignmentHorizon;
	private int recordHorizon;
	private long counter = 0;

	private static Comparator<SAMSpot> spotAlignmentStartComparator = new Comparator<SAMSpot>() {

		@Override
		public int compare(SAMSpot o1, SAMSpot o2) {
			return o1.getAlignmentStart() - o2.getAlignmentStart();
		}
	};

	private Map<String, SAMSpot> unassembledSpotsNameMap = new TreeMap<String, SAMSpot>();
	private Queue<SAMSpot> unassembledSpotsAlignmentQueue;

	private int biggestAlignmentStart = -1;

	public TemplateAssembler() {
		this(1000, 1000);
	}

	public TemplateAssembler(int alignmentHorizon, int recordHorizon) {
		this.alignmentHorizon = alignmentHorizon;
		this.recordHorizon = recordHorizon;

		unassembledSpotsAlignmentQueue = new PriorityQueue<SAMSpot>(recordHorizon, spotAlignmentStartComparator);
	}

	private void addSpot(SAMSpot spot) {
		unassembledSpotsNameMap.put(spot.getName(), spot);
		unassembledSpotsAlignmentQueue.add(spot);
	}

	private void addSingle(SAMRecord record) {
		SAMSpot spot = new SingleSAMSpot(new SAMRecordHolder(record, counter));
		addSpot(spot);
	}

	private void addPaired(SAMRecord record) {
		SAMSpot spot = new PairedSAMSpot(new SAMRecordHolder(record, counter));
		addSpot(spot);
	}

	public void addSAMRecord(SAMRecord samRecord) {
		counter++;
		biggestAlignmentStart = Math.max(biggestAlignmentStart, samRecord.getAlignmentStart());

		if (!samRecord.getProperPairFlag()) {
			addSingle(samRecord);
			return;
		}
		if (Math.abs(samRecord.getAlignmentStart() - samRecord.getMateAlignmentStart()) > alignmentHorizon) {
			addSingle(samRecord);
			return;
		}
		if (samRecord.getReadUnmappedFlag()) {
			addSingle(samRecord);
			return;
		}

		if (!unassembledSpotsNameMap.containsKey(samRecord.getReadName()))
			addPaired(samRecord);
		else {
			SAMSpot spot = unassembledSpotsNameMap.get(samRecord.getReadName());
			spot.addRecord(new SAMRecordHolder(samRecord, counter));
		}
	}

	public SAMSpot getNextAssembledTemplate() {
		SAMSpot spot = unassembledSpotsAlignmentQueue.peek();
		if (spot == null)
			return null;

		if (spot.isComplete())
			return getNextTemplate();

		if (spot.getAlignmentStart() + alignmentHorizon < biggestAlignmentStart)
			return getNextTemplate();

		if (unassembledSpotsAlignmentQueue.size() > recordHorizon)
			return getNextTemplate();

		return null;
	}

	public SAMSpot getNextTemplate() {
		SAMSpot spot = unassembledSpotsAlignmentQueue.poll();
		if (spot == null)
			return null;

		unassembledSpotsNameMap.remove(spot.getName());
		return spot;
	}
}
