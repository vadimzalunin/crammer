package uk.ac.ebi.ena.sra.cram.io;

import java.io.IOException;
import java.io.OutputStream;
import java.io.StringWriter;

public class DebuggingBitOuputStream implements BitOutputStream {
	private OutputStream os ;
	private char opSeparator = '\n';

	public DebuggingBitOuputStream(OutputStream os, char opSeparator) {
		this.os = os ;
		this.opSeparator = opSeparator;
	}

	public static void main(String[] args) throws IOException {
		StringWriter writer = new StringWriter();
		DebuggingBitOuputStream bsos = new DebuggingBitOuputStream(System.out, '\n');
		bsos.write(false);
		bsos.write(true);
		bsos.write(1, 1);
		bsos.write(10, 8);

		System.out.println(writer.toString());
	}

	@Override
	public void write(int b, int nbits) throws IOException {
		for (int i = 0; i < nbits; i++)
			append(((b >> i) & 1) == 1);
		os.write(opSeparator);
	}

	@Override
	public void write(long b, int nbits) throws IOException {
		for (int i = 0; i < nbits; i++)
			append(((b >> i) & 1) == 1);
		os.write(opSeparator);
	}

	@Override
	public void write(byte b, int nbits) throws IOException {
		for (int i = 0; i < nbits; i++)
			append(((b >> i) & 1) == 1);
		os.write(opSeparator);
	}

	@Override
	public void write(boolean bit) throws IOException {
		append(bit);
		os.write(opSeparator);
	}

	@Override
	public void write(boolean bit, long repeat) throws IOException {
		for (int i = 0; i < repeat; i++)
			append(bit);
		os.write(opSeparator);
	}

	@Override
	public void flush() throws IOException {
		os.flush();
	}

	@Override
	public void close() throws IOException {
		os.close();
	}

	private void append(boolean bit) throws IOException {
		os.write(bit ? '1' : '0');
	}
}
