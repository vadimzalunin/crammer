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

import uk.ac.ebi.ena.sra.cram.format.compression.EncodingAlgorithm;

public class Encoding implements Serializable {
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
			return foe.getAlgorithm() == getAlgorithm() && foe.getParameters().equals(getParameters());
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
