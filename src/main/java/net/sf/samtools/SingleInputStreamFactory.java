package net.sf.samtools;

import java.io.FileNotFoundException;
import java.io.InputStream;

public class SingleInputStreamFactory implements InputStreamFactory {
	private InputStream delegate;

	public SingleInputStreamFactory(InputStream delegate) {
		this.delegate = delegate;
	}

	@Override
	public InputStream createInputStream() throws FileNotFoundException,
			ExhaustedFactoryException {
		if (delegate == null)
			throw new ExhaustedFactoryException();
		InputStream is = delegate;
		delegate = null;
		return is;
	}

}
