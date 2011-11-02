package uk.ac.ebi.ena.sra.cram.index;

public class BitPointer {
	private long byteOffset;
	private byte bitOffset;

	public long getByteOffset() {
		return byteOffset;
	}

	public void setByteOffset(long byteOffset) {
		this.byteOffset = byteOffset;
	}

	public byte getBitOffset() {
		return bitOffset;
	}

	public void setBitOffset(byte bitOffset) {
		this.bitOffset = bitOffset;
	}

	@Override
	public String toString() {
		return byteOffset + "." + bitOffset;
	}
}
