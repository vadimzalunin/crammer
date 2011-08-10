package uk.ac.ebi.ena.sra.cram.format;

public class ReadAnnotation implements Comparable<ReadAnnotation> {
	// non-null
	private String key;

	public ReadAnnotation(String key) {
		if (key == null)
			throw new NullPointerException("Annotation key cannot be null.");
		this.key = key;
	}

	public String getKey() {
		return key;
	}

	@Override
	public int compareTo(ReadAnnotation o) {
		return key.compareTo(o.key);
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof ReadAnnotation))
			return false;

		return key.equals(((ReadAnnotation) obj).key);
	}
	
	@Override
	public int hashCode() {
		return key.hashCode();
	}
}
