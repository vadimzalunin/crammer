package uk.ac.ebi.ena.sra.cram.format;

public class CramReadGroup {

	private String id;
	private String sample;
	
	public CramReadGroup(String id) {
		this.id = id;
	}

	public CramReadGroup(String id, String sample) {
		this.id = id;
		this.sample = sample;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getSample() {
		return sample;
	}

	public void setSample(String sample) {
		this.sample = sample;
	}
}
