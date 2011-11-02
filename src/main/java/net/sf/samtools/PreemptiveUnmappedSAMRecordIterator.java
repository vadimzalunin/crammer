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
