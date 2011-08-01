package uk.ac.ebi.ena.sra.cram.mask;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;

public class SingleLineMaskReader implements Closeable, ReadMaskReader {
	private final BufferedReader reader;
	private final ReadMaskFactory<String> readMaskFactory;

	public SingleLineMaskReader(BufferedReader reader,
			ReadMaskFactory<String> readMaskFactory) {
		this.reader = reader;
		this.readMaskFactory = readMaskFactory;
	}

	@Override
	public PositionMask readNextMask() throws IOException,
			ReadMaskFormatException {
		String line = reader.readLine();
		if (line == null)
			return null;

		return readMaskFactory.createMask(line);
	}

	@Override
	public void close() throws IOException {
		reader.close();
	}

	// private String line;
	// private void readLine() {
	// try {
	// line = reader.readLine();
	// } catch (Exception e) {
	// throw new RuntimeException(e);
	// }
	// }
	//
	// @Override
	// public boolean hasNext() {
	// if (line != null)
	// return true;
	//
	// readLine();
	// return line != null;
	// }
	//
	// @Override
	// public PositionMask next() throws ReadMaskFormatException {
	// if (line == null)
	// readLine();
	//
	// if (line == null)
	// throw new RuntimeException("Iterator is empty.");
	//
	// PositionMask mask = readMaskFactory.createMask(line);
	// line = null;
	// return mask;
	// }

}
