package uk.ac.ebi.ena.sra.cram.impl;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;

public class ByteArraySequenceBaseProvider implements SequenceBaseProvider {
	private byte[] sequence;
	private boolean circular = false;
	private int N_extension = 0;

	public ByteArraySequenceBaseProvider(byte[] sequence) {
		this.sequence = sequence;
	}

	@Override
	public byte getBaseAt(String seqName, long position) {
		if (position < sequence.length)
			return sequence[(int) position];

		if (circular)
			return sequence[(int) (position % sequence.length)];

		if (position < sequence.length + N_extension)
			return 'N';

		throw new RuntimeException(String.format(
				"Reference position out of range: in sequence %s, length %d, position %d.", seqName, sequence.length,
				position));
	}

	@Override
	public void copyBases(String sequenceName, long from, int len, byte[] dest) throws IOException {
		try {
			if (from + len > sequence.length)
				len = sequence.length - (int) from;
			System.arraycopy(sequence, (int) from, dest, 0, len);
		} catch (ArrayIndexOutOfBoundsException e) {
			System.err.printf("Sequence: %s; from=%d; len=%d\n", sequenceName, from, len);
			throw e;
		}
	}

	public boolean isCircular() {
		return circular;
	}

	public void setCircular(boolean circular) {
		this.circular = circular;
	}

	public int getN_extension() {
		return N_extension;
	}

	public void setN_extension(int n_extension) {
		N_extension = n_extension;
	}

}
