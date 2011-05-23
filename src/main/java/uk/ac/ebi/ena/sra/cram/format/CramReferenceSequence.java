package uk.ac.ebi.ena.sra.cram.format;

public class CramReferenceSequence {
	private String name;
	private int length;

	public CramReferenceSequence() {
	}

	public CramReferenceSequence(String name, int length) {
		this.name = name;
		this.length = length;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public int getLength() {
		return length;
	}

	public void setLength(int length) {
		this.length = length;
	}

	@Override
	public String toString() {
		return "RefSeq: [name=" + name + ", length=" + length + "]";
	}
}
