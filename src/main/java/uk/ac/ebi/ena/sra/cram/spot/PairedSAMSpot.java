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

import java.util.Arrays;
import java.util.Iterator;

class PairedSAMSpot implements SAMSpot {
	private SAMRecordHolder[] records = new SAMRecordHolder[2];

	public PairedSAMSpot() {
	}

	public PairedSAMSpot(SAMRecordHolder record) {
		if (record == null)
			throw new NullPointerException("Expecting a non-null record.");

		addRecord(record);
	}

	public PairedSAMSpot(SAMRecordHolder[] records) {
		if (records == null)
			throw new NullPointerException(
					"Expecting a non-null records array.");

		if (records.length != 2)
			throw new IllegalArgumentException(
					"Expecting exactly two records. ");

		for (SAMRecordHolder record : records)
			addRecord(record);
	}

	public boolean isComplete() {
		return records[0] != null && records[1] != null;
	}

	public SAMRecordHolder getFirstMate() {
		return records[0];
	}

	public SAMRecordHolder getSecondMate() {
		return records[1];
	}

	@Override
	public Iterator<SAMRecordHolder> iterator() {
		return Arrays.asList(records).iterator();
	}

	public void setFirstMate(SAMRecordHolder record) {
		if (record == null)
			throw new NullPointerException("Expecting a non-null record.");

		if (records[0] != null)
			throw new IllegalArgumentException(
					"First of pair mates has been already set.");
		if (!record.getSamRecord().getProperPairFlag())
			throw new IllegalArgumentException("Not a proper paired record.");

		if (!record.getSamRecord().getFirstOfPairFlag()
				|| record.getSamRecord().getSecondOfPairFlag())
			throw new IllegalArgumentException("First mate of pair expected.");

		if (records[1] != null) {
			if (!record.getSamRecord().getReadName()
					.equals(records[1].getSamRecord().getReadName()))
				throw new IllegalArgumentException(
						"Read names of pair mates dont match.");
		}
		records[0] = record;
	}

	public void setSecondMate(SAMRecordHolder record) {
		if (record == null)
			throw new NullPointerException("Expecting a non-null record.");

		if (records[1] != null)
			throw new IllegalArgumentException(
					"Second of pair mates has been already set.");
		if (!record.getSamRecord().getProperPairFlag())
			throw new IllegalArgumentException("Not a proper paired record.");

		if (record.getSamRecord().getFirstOfPairFlag()
				|| !record.getSamRecord().getSecondOfPairFlag())
			throw new IllegalArgumentException("Secod mate of pair expected.");

		if (records[0] != null) {
			if (!record.getSamRecord().getReadName()
					.equals(records[0].getSamRecord().getReadName()))
				throw new IllegalArgumentException(
						"Read names of pair mates dont match.");
		}

		records[1] = record;
	}

	@Override
	public void addRecord(SAMRecordHolder record) {
		if (record.getSamRecord().getFirstOfPairFlag())
			setFirstMate(record);
		else if (record.getSamRecord().getSecondOfPairFlag())
			setSecondMate(record);
		else
			throw new IllegalArgumentException(
					"Neither first nor second of mate pairs flag is set.");
	}

	@Override
	public String getName() {
		if (records[0] != null)
			return records[0].getSamRecord().getReadName();

		if (records[1] != null)
			return records[1].getSamRecord().getReadName();

		return null;
	}

	@Override
	public int getAlignmentStart() {

		if (records[0] != null)
			return Math.min(records[0].getSamRecord().getAlignmentStart(),
					records[0].getSamRecord().getMateAlignmentStart());

		if (records[1] != null)
			return Math.min(records[1].getSamRecord().getAlignmentStart(),
					records[1].getSamRecord().getMateAlignmentStart());

		return -1;
	}

}
