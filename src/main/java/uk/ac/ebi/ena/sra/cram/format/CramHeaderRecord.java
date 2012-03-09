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
