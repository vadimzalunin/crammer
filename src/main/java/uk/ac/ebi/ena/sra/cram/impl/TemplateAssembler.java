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

import java.util.Map;
import java.util.Queue;

import net.sf.samtools.SAMRecord;

public class TemplateAssembler {
	public static final int DEFAULT_ALIGNMENT_HORIZON = 10000;
	public static final int DEFAULT_RECORD_HORIZON = 10000;

	private static long recordCounter = 0;

	private Queue<Record> alignmentQueue;
	private Map<String, Template> templatesByName;

	private int alignmentHorizon;
	private int recordHorizon;

	private int lastAddedAlignmentStart = 0;

	public TemplateAssembler() {
		this(DEFAULT_ALIGNMENT_HORIZON, DEFAULT_RECORD_HORIZON);
	}

	public TemplateAssembler(int alignmentHorizon, int recordHorizon) {
		this.alignmentHorizon = alignmentHorizon;
		this.recordHorizon = recordHorizon;
	}

	/**
	 * Expecting valid input: alignment order, only one read for single-layout
	 * template. Single-layout templates are not cached.
	 * 
	 * @param samRecord
	 */
	public void addSAMRecord(SAMRecord samRecord) {

		// autofix unmapped out-of-alignment reads:
		if (samRecord.getReadUnmappedFlag()
				&& samRecord.getAlignmentStart() < lastAddedAlignmentStart)
			samRecord.setAlignmentStart(lastAddedAlignmentStart);

		Template template = templatesByName.get(samRecord.getReadName());
		if (template == null) {
			template = new Template();
			Record record = new Record(samRecord, template);

			if (!samRecord.getReadPairedFlag()) {
				template.firstRecord = record;
				template.lastRecord = record;
			} else {
				if (samRecord.getFirstOfPairFlag())
					template.firstRecord = record;
				else
					template.lastRecord = record;

				if (template.isComplete())
					templatesByName.remove(samRecord.getReadName());
				else
					templatesByName.put(samRecord.getReadName(), template);
			}
			alignmentQueue.add(record);
			if (lastAddedAlignmentStart < samRecord.getAlignmentStart())
				lastAddedAlignmentStart = samRecord.getAlignmentStart();

			return;
		}

		Record record = new Record(samRecord, template);
		if (samRecord.getFirstOfPairFlag())
			template.firstRecord = record;
		else
			template.lastRecord = record;

		alignmentQueue.add(record);
	}

	public Record poll() {
		Record nextRecordOnAlignment = alignmentQueue.peek();
		if (nextRecordOnAlignment == null)
			return null;

		if (nextRecordOnAlignment.template == null
				|| nextRecordOnAlignment.template.isSingle()
				|| nextRecordOnAlignment.template.isComplete())
			return alignmentQueue.poll();

		if ((lastAddedAlignmentStart
				- nextRecordOnAlignment.samRecord.getAlignmentStart() > alignmentHorizon)
				|| (alignmentQueue.size() > recordHorizon))
			return takeNextRecordWithIncompleteTemplate();

		return null;
	}

	private Record takeNextRecordWithIncompleteTemplate() {
		Record nextRecordOnAlignment = alignmentQueue.poll();
		divorceSAMRecord(nextRecordOnAlignment);
		templatesByName.remove(nextRecordOnAlignment.samRecord.getReadName());
		return nextRecordOnAlignment;
	}

	private void divorceSAMRecord(Record record) {
		record.samRecord.setReadPairedFlag(false);
		record.samRecord
				.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
		record.samRecord.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
		record.samRecord.setInferredInsertSize(0);

		record.template = null;
		record.nextRecordInTemplate = null;
		record.prevRecordInTemplate = null;
	}

	public class Template {
		private Record firstRecord;
		private Record lastRecord;

		public boolean isSingle() {
			return firstRecord != null && firstRecord == lastRecord;
		}

		public boolean isComplete() {
			Record record = firstRecord;
			while (record != null && record != lastRecord) {
				record = record.nextRecordInTemplate;
			}
			return record != null;
		}

		// public long getMinAlignmentStart() {
		// long min = Long.MAX_VALUE;
		// Record record = firstRecord;
		// while (record != null && record != lastRecord) {
		// min = Math.min(min, record.samRecord.getAlignmentStart());
		// record = record.nextRecordInTemplate;
		// }
		//
		// return min;
		// }
		//
		// public long getMaxAlignmentStart() {
		// long max = Long.MIN_VALUE;
		// Record record = firstRecord;
		// while (record != null && record != lastRecord) {
		// max = Math.max(max, record.samRecord.getAlignmentStart());
		// record = record.nextRecordInTemplate;
		// }
		//
		// return max;
		// }
	}

	public class Record {
		private long inStreamRecordPosition;
		private Record prevRecordInTemplate;
		private Record nextRecordInTemplate;
		private Template template;

		private SAMRecord samRecord;

		public Record(SAMRecord samRecord, Template template) {
			this.samRecord = samRecord;
			this.template = template;
			inStreamRecordPosition = recordCounter++;
		}

		public long getDistanceToNextFragment() {
			if (nextRecordInTemplate == null)
				return -1;
			return nextRecordInTemplate.inStreamRecordPosition
					- inStreamRecordPosition;
		}

		public String getName() {
			return samRecord.getReadName();
		}

		public long getAlignmentStart() {
			return samRecord.getAlignmentStart();
		}

		public Record getPrevRecordInTemplate() {
			return prevRecordInTemplate;
		}

		public Record getNextRecordInTemplate() {
			return nextRecordInTemplate;
		}
	}
}
