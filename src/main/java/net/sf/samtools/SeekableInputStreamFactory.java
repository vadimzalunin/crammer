package net.sf.samtools;

import java.io.IOException;

import net.sf.samtools.util.SeekableStream;

public interface SeekableInputStreamFactory {

	public SeekableStream createSeekableStream() throws IOException,
			ExhaustedFactoryException;
}
