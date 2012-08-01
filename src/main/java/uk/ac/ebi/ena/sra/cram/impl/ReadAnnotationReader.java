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
package uk.ac.ebi.ena.sra.cram.impl;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import uk.ac.ebi.ena.sra.cram.CramException;
import uk.ac.ebi.ena.sra.cram.format.ReadAnnotation;

public class ReadAnnotationReader {

	private static final String DLM = "\t";
	private Set<ReadAnnotation> annnotations;
	private BufferedReader bufferedReader;
	private Map<String, ReadAnnotation> key2ReadAnnotationMap = new HashMap<String, ReadAnnotation>();

	public ReadAnnotationReader(BufferedReader bufferedReader)
			throws IOException {
		this.bufferedReader = bufferedReader;
		String dictionaryLine = bufferedReader.readLine();
		for (String key : dictionaryLine.split(DLM))
			key2ReadAnnotationMap.put(key, new ReadAnnotation(key));

		this.annnotations = Collections
				.unmodifiableSet(new TreeSet<ReadAnnotation>(
						key2ReadAnnotationMap.values()));
	}

	public List<ReadAnnotation> listUniqAnnotations() {
		return Collections.unmodifiableList(new ArrayList<ReadAnnotation>(
				key2ReadAnnotationMap.values()));
	}

	public Collection<ReadAnnotation> nextReadAnnotations() throws IOException,
			ReadAnnotationException {
		String line = bufferedReader.readLine();
		if (line == null)
			return null;

		List<ReadAnnotation> annoList = new ArrayList<ReadAnnotation>();
		for (String chunk : line.split(DLM))
			if (chunk.length() > 0)
				annoList.add(createReadAnnotation(chunk));

		return annoList;
	}

	private ReadAnnotation createReadAnnotation(String chunk)
			throws ReadAnnotationException {
		ReadAnnotation readAnnotation = key2ReadAnnotationMap.get(chunk);
		if (readAnnotation == null)
			throw new ReadAnnotationException(
					"Annotaion not found in the dictionary: " + chunk);

		return readAnnotation;
	}

	public static class ReadAnnotationException extends CramException {

		public ReadAnnotationException(String message) {
			super(message);
		}

	}
}
