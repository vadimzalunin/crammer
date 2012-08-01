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
package uk.ac.ebi.ena.sra.cram.format;

import java.io.Serializable;
import java.util.Collection;

public class CramRecordBlock implements Serializable {
	private transient Collection<CramRecord> records;
	private String sequenceName;
	private int sequenceLength;
	private long firstRecordPosition;
	private long recordCount;
	private int readLength;
	private CramCompression compression;
	private boolean positiveStrandBasePositionReversed = false;
	private boolean negativeStrandBasePositionReversed = !positiveStrandBasePositionReversed;

	private boolean unmappedReadQualityScoresIncluded = false;
	private boolean substitutionQualityScoresIncluded = true;
	private boolean maskedQualityScoresIncluded = true;
	
	public boolean losslessQualityScores = false;
	public boolean preserveReadNames = false ;
	
	public String getSequenceName() {
		return sequenceName;
	}

	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}

	public Collection<CramRecord> getRecords() {
		return records;
	}

	public void setRecords(Collection<CramRecord> records) {
		this.records = records;
	}

	public long getFirstRecordPosition() {
		return firstRecordPosition;
	}

	public void setFirstRecordPosition(long firstRecordPosition) {
		this.firstRecordPosition = firstRecordPosition;
	}

	public long getRecordCount() {
		return recordCount;
	}

	public void setRecordCount(long recordCount) {
		this.recordCount = recordCount;
	}

	public int getReadLength() {
		return readLength;
	}

	@Deprecated
	public void setReadLength(int readLength) {
		this.readLength = readLength;
	}

	public CramCompression getCompression() {
		return compression;
	}

	public void setCompression(CramCompression compression) {
		this.compression = compression;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof CramRecordBlock))
			return false;

		CramRecordBlock foe = (CramRecordBlock) obj;

		if (sequenceName == null)
			if (foe.sequenceName != null)
				return false;
		if (!sequenceName.equals(foe.sequenceName))
			return false;
		if (sequenceLength != foe.sequenceLength)
			return false;
		if (firstRecordPosition != foe.firstRecordPosition)
			return false;
		if (recordCount != foe.recordCount)
			return false;
		if (recordCount != foe.recordCount)
			return false;
		if (readLength != foe.readLength)
			return false;
		if (positiveStrandBasePositionReversed != foe.positiveStrandBasePositionReversed)
			return false;
		if (negativeStrandBasePositionReversed != foe.negativeStrandBasePositionReversed)
			return false;
		if (unmappedReadQualityScoresIncluded != foe.unmappedReadQualityScoresIncluded)
			return false;
		if (substitutionQualityScoresIncluded != foe.substitutionQualityScoresIncluded)
			return false;
		if (maskedQualityScoresIncluded != foe.maskedQualityScoresIncluded)
			return false;

		if (!compression.equals(foe.compression))
			return false;

		return true;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(getClass().getSimpleName()).append(" [\n");
		sb.append("sequenceName=").append(sequenceName).append("\n");
		sb.append("sequenceLength=").append(sequenceLength).append("\n");
		sb.append("firstRecordPosition=").append(firstRecordPosition)
				.append("\n");
		sb.append("recordCount=").append(recordCount).append("\n");
		sb.append("readLength=").append(readLength).append("\n");
		sb.append("unmappedReadQualityScoresIncluded=")
				.append(unmappedReadQualityScoresIncluded).append("\n");
		sb.append("substitutionQualityScoresIncluded=")
				.append(substitutionQualityScoresIncluded).append("\n");
		sb.append("maskedQualityScoresIncluded=")
				.append(maskedQualityScoresIncluded).append("\n");
		sb.append("compression=").append(compression).append("]");
		return sb.toString();
	}

	public boolean isPositiveStrandBasePositionReversed() {
		return positiveStrandBasePositionReversed;
	}

	public void setPositiveStrandBasePositionReversed(
			boolean positiveStrandBasePositionReversed) {
		this.positiveStrandBasePositionReversed = positiveStrandBasePositionReversed;
	}

	public boolean isNegativeStrandBasePositionReversed() {
		return negativeStrandBasePositionReversed;
	}

	public void setNegativeStrandBasePositionReversed(
			boolean negativeStrandBasePositionReversed) {
		this.negativeStrandBasePositionReversed = negativeStrandBasePositionReversed;
	}

	public int getSequenceLength() {
		return sequenceLength;
	}

	public void setSequenceLength(int sequenceLength) {
		this.sequenceLength = sequenceLength;
	}

	public boolean isUnmappedReadQualityScoresIncluded() {
		return unmappedReadQualityScoresIncluded;
	}

	public void setUnmappedReadQualityScoresIncluded(
			boolean unmappedReadQualityScoresIncluded) {
		this.unmappedReadQualityScoresIncluded = unmappedReadQualityScoresIncluded;
	}

	public boolean isSubstitutionQualityScoresIncluded() {
		return substitutionQualityScoresIncluded;
	}

	public void setSubstitutionQualityScoresIncluded(
			boolean substitutionQualityScoresIncluded) {
		this.substitutionQualityScoresIncluded = substitutionQualityScoresIncluded;
	}

	public boolean isMaskedQualityScoresIncluded() {
		return maskedQualityScoresIncluded;
	}

	public void setMaskedQualityScoresIncluded(
			boolean maskedQualityScoresIncluded) {
		this.maskedQualityScoresIncluded = maskedQualityScoresIncluded;
	}

}
