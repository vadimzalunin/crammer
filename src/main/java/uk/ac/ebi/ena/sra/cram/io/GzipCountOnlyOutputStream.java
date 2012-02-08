package uk.ac.ebi.ena.sra.cram.io;

import java.io.IOException;
import java.io.OutputStream;
import java.util.zip.Deflater;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.compress.compressors.bzip2.BZip2CompressorOutputStream;

public class GzipCountOnlyOutputStream extends OutputStream {
	public enum COMPRESSOR {
		NONE, GZIP, BZIP2, LZMA;
	}

	private CountOnlyOutputStream countOS = new CountOnlyOutputStream();
	private OutputStream gzipOS;

	public GzipCountOnlyOutputStream(COMPRESSOR compressor) {
		try {
			switch (compressor) {
			case NONE:
				gzipOS = countOS;
				break;
			case GZIP:
				gzipOS = new GZIPOutputStream(countOS) {
					{
//						def.setLevel(Deflater.BEST_COMPRESSION);
					}
				};
				break;
			case BZIP2:
				gzipOS = new BZip2CompressorOutputStream(countOS);
				break;
			default:
				throw new RuntimeException("Unsupported compression: " + compressor);
			}

		} catch (IOException e) {
			throw new RuntimeException("Failed to create gzip output stream.", e);
		}
	}

	public GzipCountOnlyOutputStream() {
		this (COMPRESSOR.GZIP) ;
	}

	public final void writeInt(int v) throws IOException {
		write((v >>> 24) & 0xFF);
		write((v >>> 16) & 0xFF);
		write((v >>> 8) & 0xFF);
		write((v >>> 0) & 0xFF);
	}

	public final void writeFloat(float v) throws IOException {
		writeInt(Float.floatToIntBits(v));
	}

	public final void writeByte(int v) throws IOException {
		write(v);
	}

	public final void writeShort(int v) throws IOException {
		write((v >>> 8) & 0xFF);
		write((v >>> 0) & 0xFF);
	}

	@Override
	public void write(int b) throws IOException {
		gzipOS.write(b);
	}

	@Override
	public void write(byte[] b) throws IOException {
		gzipOS.write(b);
	}

	@Override
	public void write(byte[] b, int off, int len) throws IOException {
		gzipOS.write(b, off, len);
	}

	@Override
	public void close() throws IOException {
		gzipOS.close();
	}

	@Override
	public void flush() throws IOException {
		gzipOS.flush();
	}

	public long getCount() {
		return countOS.getCount();
	}

}
