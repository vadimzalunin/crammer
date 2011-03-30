package uk.ac.ebi.ena.sra.cram.format;

import java.io.Serializable;
import java.util.Collection;

public class CramRecordBlock implements Serializable {
	String sequenceName;
	transient Collection<CramRecord> records;
	int firstRecordPosition;
	int recordCount;
	int readLength;
	CramCompression compression ;

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

	public int getFirstRecordPosition() {
		return firstRecordPosition;
	}

	public void setFirstRecordPosition(int firstRecordPosition) {
		this.firstRecordPosition = firstRecordPosition;
	}

	public int getRecordCount() {
		return recordCount;
	}

	public void setRecordCount(int recordCount) {
		this.recordCount = recordCount;
	}

	public int getReadLength() {
		return readLength;
	}

	public void setReadLength(int readLength) {
		this.readLength = readLength;
	}

	public CramCompression getCompression() {
		return compression;
	}

	public void setCompression(CramCompression compression) {
		this.compression = compression;
	}
}
