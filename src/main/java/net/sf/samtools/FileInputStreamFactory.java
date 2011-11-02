package net.sf.samtools;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;

public class FileInputStreamFactory implements InputStreamFactory {
	private File file;
	private boolean buffered = true;

	public FileInputStreamFactory(File file) {
		this.file = file;
	}

	@Override
	public InputStream createInputStream() throws FileNotFoundException {
		if (buffered)
			return new BufferedInputStream(new FileInputStream(file));
		else
			return new FileInputStream(file);
	}

}
