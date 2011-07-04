package uk.ac.ebi.ena.sra.cram.format;

import java.io.Serializable;
import java.util.Collection;
import java.util.List;

public class CramHeader implements Serializable {
	private String version;
	private Collection<CramReferenceSequence> referenceSequences;

	public String getVersion() {
		return version;
	}

	public void setVersion(String version) {
		this.version = version;
	}

	public Collection<CramReferenceSequence> getReferenceSequences() {
		return referenceSequences;
	}

	public void setReferenceSequences(
			Collection<CramReferenceSequence> referenceSequences) {
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
}
