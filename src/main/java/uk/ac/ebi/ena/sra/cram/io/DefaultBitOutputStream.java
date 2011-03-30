package uk.ac.ebi.ena.sra.cram.io;

import java.io.IOException;
import java.io.OutputStream;

public class DefaultBitOutputStream extends OutputStream implements
		BitOutputStream {

	private static final byte[] bitMasks = new byte[8];
	static {
		for (byte i = 0; i < 8; i++)
			bitMasks[i] = (byte) (~(0xFF >>> i));
	}

	private final OutputStream out;

	private int bufferByte = 0;
	private int bufferedNumberOfBits = 0;

	public DefaultBitOutputStream(OutputStream delegate) {
		this.out = delegate;
	}

	@Override
	public void write(int value) throws IOException {
		write(toBytes(value));
	}

	private final static byte[] toBytes(int value) {
		final byte[] bytes = new byte[4];
		bytes[0] = (byte) (value >>> 24);
		bytes[1] = (byte) (value >>> 16);
		bytes[2] = (byte) (value >>> 8);
		bytes[3] = (byte) (value >>> 0);
		return bytes;
	}

	public void write(int value, int nofBitsToWrite) throws IOException {
		if (nofBitsToWrite < 1 || nofBitsToWrite > 32)
			throw new IOException("Expecting 1 to 32 bits.");

		if (nofBitsToWrite < 8)
			writeBits(value, nofBitsToWrite);
		else {
			for (int i = 0;; i += 8) {
				final int v = value >>> (24 - i);
				if (i >= nofBitsToWrite) {
					writeBits(v, i % 8);
					break;
				} else
					writeBits(v, 8);
			}
		}
	}

	public void writeByte(int value) throws IOException {
		if (bufferedNumberOfBits == 0)
			putByte(value);
		else {
			bufferByte = (value & 0xFF >>> bufferedNumberOfBits) | bufferByte;
			putByte(bufferByte);
			bufferByte = (value << (8 - bufferedNumberOfBits)) & 0xFF;
		}
	}

	public void writeBits(int value, int nofBitsToWrite) throws IOException {
		if (nofBitsToWrite < 0 || nofBitsToWrite > 8) 
			throw new IOException("Expecting 0 to 8 bits.");
		
		if (nofBitsToWrite == 8)
			writeByte(value);
		else {
			if (bufferedNumberOfBits == 0) {
				bufferByte = (value << (8 - nofBitsToWrite)) & 0xFF;
				bufferedNumberOfBits = nofBitsToWrite;
			} else {
				value = (byte) (value & ~bitMasks[8 - nofBitsToWrite]);
				int bits = 8 - bufferedNumberOfBits - nofBitsToWrite;
				if (bits < 0) {
					bits = -bits;
					bufferByte |= (value >>> bits);
					putByte(bufferByte);
					bufferByte = (value << (8 - bits)) & 0xFF;
					bufferedNumberOfBits = bits;
				} else if (bits == 0) {
					bufferByte = bufferByte | value;
					putByte(bufferByte);
					bufferedNumberOfBits = 0;
				} else {
					bufferByte = bufferByte | (value << bits);
					bufferedNumberOfBits = 8 - bits;
				}
			}
		}
	}

	@Override
	public void close() throws IOException {
		flush();
		out.close();
	}

	@Override
	public void flush() throws IOException {
		if (bufferedNumberOfBits > 0)
			putByte(bufferByte);

		out.flush();
	}

	private void putByte(int value) throws IOException {
		out.write(value);
	}
}
