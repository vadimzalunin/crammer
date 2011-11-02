package net.sf.samtools;

import java.io.IOException;

import net.sf.samtools.util.SeekableStream;

public class SingleSeekableStreamFactory implements SeekableInputStreamFactory {
	private SeekableStream delegate;

	public SingleSeekableStreamFactory(SeekableStream delegate) {
		this.delegate = delegate;
	}

	@Override
	public SeekableStream createSeekableStream() throws IOException,
			ExhaustedFactoryException {
		SeekableStream ss = delegate;
		delegate = null;
		return ss;
	}

}
