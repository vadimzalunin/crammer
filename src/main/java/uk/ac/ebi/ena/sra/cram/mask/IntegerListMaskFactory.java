package uk.ac.ebi.ena.sra.cram.mask;

import java.nio.IntBuffer;

public class IntegerListMaskFactory implements ReadMaskFactory<String> {
	public static final String DEFAULT_DEMLIITER = " ";
	public static final int DEFAULT_BUFFER_SIZE = 1024;

	private String delimiter;
	private IntBuffer buf;

	public IntegerListMaskFactory(String delimiter, int bufSize) {
		if (delimiter == null)
			throw new NullPointerException("Delimiter is null.");
		if (bufSize < 1)
			throw new IllegalAccessError("Buffer size must be greater then 0.");

		this.delimiter = delimiter;
		this.buf = IntBuffer.allocate(bufSize);
	}

	public IntegerListMaskFactory(String delimiter) {
		this(delimiter, DEFAULT_BUFFER_SIZE);
	}

	public IntegerListMaskFactory() {
		this(DEFAULT_DEMLIITER, DEFAULT_BUFFER_SIZE);
	}

	@Override
	public PositionMask createMask(String line) throws ReadMaskFormatException {
		if (line.length() == 0)
			return ArrayPositionMask.EMPTY_INSTANCE;
		
		buf.clear() ;
		try {
			for (String chunk : line.split(delimiter))
				buf.put(Integer.valueOf(chunk));
		} catch (NumberFormatException e) {
			throw new ReadMaskFormatException("Expecting integers in "
					+ line.substring(0, Math.min(10, line.length())), e);
		}
		buf.flip();
		int[] array = new int[buf.limit()] ;
		buf.get(array) ;
		PositionMask mask = new ArrayPositionMask(array);

		return mask;
	}

}
