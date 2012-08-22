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
package uk.ac.ebi.ena.sra.cram.index;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class AlignmentIndex implements Iterable<RecordPointer> {

	private List<RecordPointer> pointers = new ArrayList<RecordPointer>();
	private RecordPointer searchPointer = new RecordPointer();

	public RecordPointer findRecordPointer(long alStart) {
		searchPointer.setAlignmentStart(alStart);
		int pointerIndex = Collections.binarySearch(pointers, searchPointer);

		if (pointerIndex == 0 || pointerIndex == -1)
			return pointers.get(0);

		if (pointerIndex > 0)
			return pointers.get(pointerIndex - 1);

		pointerIndex = -2 - pointerIndex;

		// if (pointerIndex < pointers.size())
		return pointers.get(pointerIndex);
		// else
		// return pointers.get(pointerIndex - 1);
	}

	void addRecordPointer(RecordPointer pointer) {
		pointers.add(pointer);
	}

	@Override
	public Iterator<RecordPointer> iterator() {
		return pointers.iterator();
	}
}
