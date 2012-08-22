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

public class ReadAnnotation implements Comparable<ReadAnnotation> {
	// non-null
	private String key;

	public ReadAnnotation(String key) {
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

		ReadAnnotation foe = (ReadAnnotation) obj;
		return key.equals(foe.key);
	}

	@Override
	public int hashCode() {
		return key.hashCode();
	}

}
