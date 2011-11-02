package net.sf.samtools;

import java.io.IOException;
import java.io.InputStream;

public interface InputStreamFactory {

	public InputStream createInputStream() throws IOException,
			ExhaustedFactoryException;
}
