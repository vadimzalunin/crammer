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

import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

public class CramHeaderRecord {
	private String tag;
	private Map<String, String> content = new TreeMap<String, String>();

	public CramHeaderRecord(String tag) {
		this.tag = tag;
	}

	public String getTag() {
		return tag;
	}

	public String getValue(String key) {
		return content.get(key);
	}

	public void setValue(String key, String value) {
		content.put(key, value);
	}

	public void setValueIfNotNull(String key, String value) {
		if (value != null)
			content.put(key, value);
	}

	public Set<String> getKeySet() {
		return content.keySet();
	}
}
