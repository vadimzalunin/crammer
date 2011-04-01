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

		if (nofBitsToWrite <= 8)
			write((byte)value, nofBitsToWrite);
		else {
			for (int i = 0;; i += 8) {
				final int v = value >>> (24 - i);
				if (i >= nofBitsToWrite) {
					write((byte)v, i % 8);
					break;
				} else
					write((byte)v, 8);
			}
		}
	}

	public void writeByte(int value) throws IOException {
		if (bufferedNumberOfBits == 0)
			out.write(value);
		else {
			bufferByte = (value & 0xFF >>> bufferedNumberOfBits) | bufferByte;
			out.write(bufferByte);
			bufferByte = (value << (8 - bufferedNumberOfBits)) & 0xFF;
		}
	}

	public void write(byte value, int nofBitsToWrite) throws IOException {
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
					out.write(bufferByte);
					bufferByte = (value << (8 - bits)) & 0xFF;
					bufferedNumberOfBits = bits;
				} else if (bits == 0) {
					bufferByte = bufferByte | value;
					out.write(bufferByte);
					bufferedNumberOfBits = 0;
				} else {
					bufferByte = bufferByte | (value << bits);
					bufferedNumberOfBits = 8 - bits;
				}
			}
		}
	}

	public void write(long value, int nofBitsToWrite) throws IOException {
		if (nofBitsToWrite < 1 || nofBitsToWrite > 64)
			throw new IOException("Expecting 1 to 64 bits.");

		if (nofBitsToWrite < 33)
			write((int) value, nofBitsToWrite);
		else {
			int highBits = (int) (value >>> 32);
			int lowBits = (int) value;

			write(highBits, nofBitsToWrite - 32);
			write(lowBits, 32);
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
			out.write(bufferByte);

		out.flush();
	}


}
