package uk.ac.ebi.ena.sra.cram.stats;

import java.util.Map;
import java.util.TreeMap;

import uk.ac.ebi.ena.sra.cram.format.ReadTag;

public class ReadTagStats {

	private Map<String, LimitedBag> tagBags = new TreeMap<String, LimitedBag>();

	public void add(ReadTag tag) {
		if (!tagBags.containsKey(tag.getKey())) 
			tagBags.put(tag.getKey(), new LimitedBag(tag.getKey()));
		
		tagBags.get(tag.getKey()).add(tag.getValue());
	}
}
