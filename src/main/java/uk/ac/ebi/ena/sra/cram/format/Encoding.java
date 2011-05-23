package uk.ac.ebi.ena.sra.cram.format;

import java.io.Serializable;

import uk.ac.ebi.ena.sra.cram.format.compression.EncodingAlgorithm;

public class Encoding implements Serializable{
	private EncodingAlgorithm algorithm;
	private String parameters;

	public Encoding() {
	}

	public Encoding(EncodingAlgorithm algorithm, String parameters) {
		super();
		this.setAlgorithm(algorithm);
		this.setParameters(parameters);
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof Encoding) {
			Encoding foe = (Encoding) obj;
			return foe.getAlgorithm() == getAlgorithm()
					&& foe.getParameters().equals(getParameters());
		} else
			return false;
	}
	
	@Override
	public String toString() {
		return getClass().getSimpleName() + " [" + getAlgorithm().name() + ", " + getParameters() + "]";
	}

	public void setAlgorithm(EncodingAlgorithm algorithm) {
		this.algorithm = algorithm;
	}

	public EncodingAlgorithm getAlgorithm() {
		return algorithm;
	}

	public void setParameters(String parameters) {
		this.parameters = parameters;
	}

	public String getParameters() {
		return parameters;
	}
}
