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
package net.sf.samtools;

import net.sf.samtools.util.CloseableIterator;

public class PreemptiveSAMRecordIterator implements
		CloseableIterator<SAMRecord> {
	CloseableIterator<SAMRecord> delegate;
	SAMRecord nextRecord;
	boolean eof = false;
	long start = 0;
	long end = Long.MAX_VALUE;
	boolean contained = false;

	public PreemptiveSAMRecordIterator(CloseableIterator<SAMRecord> delegate) {
		this(delegate, 0, Long.MAX_VALUE, false);
	}

	public PreemptiveSAMRecordIterator(CloseableIterator<SAMRecord> delegate,
			long start) {
		this (delegate, start, Long.MAX_VALUE, false); 
	}

	public PreemptiveSAMRecordIterator(CloseableIterator<SAMRecord> delegate,
			long start, long end, boolean contained) {
		this.delegate = delegate;
		this.start = start;
		this.end = end;
		this.contained = contained;
	}

	@Override
	public boolean hasNext() {
		if (eof)
			return false;

		if (nextRecord == null) {
			// gettting first record:
			advance() ;
		}
		return !eof;
	}

	@Override
	public SAMRecord next() {
		if (eof)
			throw new RuntimeException("No more records.");

		SAMRecord record = nextRecord;
		advance();
		return record;
	}

	protected SAMRecord advance() {
		if (!delegate.hasNext()) {
			eof = true;
			return nextRecord;
		}

		do {
			nextRecord = delegate.next();
		} while ((contained && nextRecord.getAlignmentStart() < start)
				|| (!contained && nextRecord.getAlignmentStart()
						+ nextRecord.getReadLength() < start));

		if ((!contained && nextRecord.getAlignmentStart() > end)
				|| (contained && nextRecord.getAlignmentEnd() > end))
			eof = true;

		return nextRecord;
	}

	@Override
	public void remove() {
		delegate.remove();
	}

	@Override
	public void close() {
		delegate.close();
	}

}
