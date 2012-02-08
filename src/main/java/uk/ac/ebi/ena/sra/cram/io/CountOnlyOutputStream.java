package uk.ac.ebi.ena.sra.cram.io;

import java.io.IOException;
import java.io.OutputStream;

public class CountOnlyOutputStream extends OutputStream {
	private long count = 0L;

	@Override
	public void write(int b) throws IOException {
		count++;
	}

	@Override
	public void write(byte[] b) throws IOException {
		count += b.length;
	}

	@Override
	public void write(byte[] b, int off, int len) throws IOException {
		count += len;
	}

	public long getCount() {
		return count;
	}

}
