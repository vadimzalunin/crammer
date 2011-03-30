package uk.ac.ebi.ena.sra.cram.format;

import java.io.Serializable;
import java.util.Collection;

public class CramHeader implements Serializable{
	String version;
	Collection<CramRecordBlock> blocks;

	public String getVersion() {
		return version;
	}

	public void setVersion(String version) {
		this.version = version;
	}

	public Collection<CramRecordBlock> getBlocks() {
		return blocks;
	}

	public void setBlocks(Collection<CramRecordBlock> blocks) {
		this.blocks = blocks;
	}

}
