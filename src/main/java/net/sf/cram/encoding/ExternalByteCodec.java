package net.sf.cram.encoding;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import net.sf.cram.io.BitInputStream;
import net.sf.cram.io.BitOutputStream;


public class ExternalByteCodec implements BitCodec<Byte> {
	private OutputStream os;
	private InputStream is;

	public ExternalByteCodec(OutputStream os, InputStream is) {
		this.os = os;
		this.is = is;
	}

	@Override
	public Byte read(BitInputStream bis) throws IOException {
		return (byte) is.read();
	}

	@Override
	public long write(BitOutputStream bos, Byte object) throws IOException {
		os.write(object);
		return 8;
	}

	@Override
	public long numberOfBits(Byte object) {
		return 8;
	}

	@Override
	public Byte read(BitInputStream bis, int len) throws IOException {
		throw new RuntimeException("Not implemented.") ;
	}
}
