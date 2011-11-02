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
