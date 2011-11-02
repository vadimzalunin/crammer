package net.sf.samtools;

import java.io.File;
import java.io.FileNotFoundException;

import net.sf.samtools.util.SeekableFileStream;
import net.sf.samtools.util.SeekableStream;

public class SeekableFileStreamFactory implements SeekableInputStreamFactory {
	private File file;

	public SeekableFileStreamFactory(File file) {
		super();
		this.file = file;
	}

	@Override
	public SeekableStream createSeekableStream() throws FileNotFoundException {
		return new SeekableFileStream(file);
	}

}
