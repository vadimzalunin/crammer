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

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.collections.map.MultiValueMap;

public class CramHeader implements Serializable {
	private String version;
	private List<CramReferenceSequence> referenceSequences;
	private List<ReadAnnotation> readAnnotations;
	private List<CramReadGroup> readGroups;

	private List<CramHeaderRecord> records = new ArrayList<CramHeaderRecord>();

	public static final String VERSION;
	static {
		if (CramHeader.class.getPackage().getImplementationVersion() != null)
			VERSION = CramHeader.class.getPackage().getImplementationVersion();
		else {
			// screw this, hard coding for now:
			VERSION = "0.9";
		}
	}

	public CramHeader() {
		if (version == null)
			version = VERSION;
	}

	public String getMajorVersion() {
		return version.substring(0, 3);
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

	public void setReferenceSequences(List<CramReferenceSequence> referenceSequences) {
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

	public void setReadAnnotations(List<ReadAnnotation> readAnnotationDictionary) {
		this.readAnnotations = readAnnotationDictionary;
	}

	public List<CramReadGroup> getReadGroups() {
		return readGroups;
	}

	public void setReadGroups(List<CramReadGroup> readGroups) {
		this.readGroups = readGroups;
	}

	public List<CramHeaderRecord> getRecords() {
		return records;
	}

	public void setRecords(List<CramHeaderRecord> records) {
		this.records = records;
	}

}
