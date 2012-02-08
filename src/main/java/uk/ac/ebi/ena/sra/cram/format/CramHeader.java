package uk.ac.ebi.ena.sra.cram.format;

import java.io.Serializable;
import java.util.List;

public class CramHeader implements Serializable {
	private String version;
	private List<CramReferenceSequence> referenceSequences;
	private List<ReadAnnotation> readAnnotations;
	private List<CramReadGroup> readGroups;
	
	public CramHeader() {
		version = "0.6" ;
	}

	public String getVersion() {
		return version;
	}

	public void setVersion(String version) {
		this.version = version;
	}

	public List<CramReferenceSequence> getReferenceSequences() {
		return referenceSequences;
	}

	public void setReferenceSequences(
			List<CramReferenceSequence> referenceSequences) {
		this.referenceSequences = referenceSequences;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("version=").append(version);
		if (referenceSequences == null || referenceSequences.isEmpty())
			return sb.toString();

		for (CramReferenceSequence s : referenceSequences) {
			sb.append(", ");
			sb.append(s.toString());
		}
		return sb.toString();
	}

	public List<ReadAnnotation> getReadAnnotations() {
		return readAnnotations;
	}

	public void setReadAnnotations(
			List<ReadAnnotation> readAnnotationDictionary) {
		this.readAnnotations = readAnnotationDictionary;
	}

	public List<CramReadGroup> getReadGroups() {
		return readGroups;
	}

	public void setReadGroups(List<CramReadGroup> readGroups) {
		this.readGroups = readGroups;
	}
}
