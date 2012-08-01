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

public class PreemptiveUnmappedSAMRecordIterator implements
		CloseableIterator<SAMRecord> {
	CloseableIterator<SAMRecord> delegate;
	SAMRecord nextRecord;
	boolean eof = false;
	
	public PreemptiveUnmappedSAMRecordIterator(
			CloseableIterator<SAMRecord> delegate) {
		super();
		this.delegate = delegate;
	}

	@Override
	public boolean hasNext() {
		if (eof)
			return false;

		if (nextRecord == null)
			advance();
		return !eof;
	}

	@Override
	public SAMRecord next() {
		if (eof)
			throw new RuntimeException("No more records.");

		return advance();
	}

	protected SAMRecord advance() {
		while (delegate.hasNext()) {
			nextRecord = delegate.next();
			if (nextRecord.getReadUnmappedFlag())
				return nextRecord;
		}

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
