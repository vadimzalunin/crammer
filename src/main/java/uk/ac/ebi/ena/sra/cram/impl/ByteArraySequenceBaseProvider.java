package uk.ac.ebi.ena.sra.cram.impl;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;

public class ByteArraySequenceBaseProvider implements SequenceBaseProvider {
	private byte[] sequence;

	public ByteArraySequenceBaseProvider(byte[] sequence) {
		super();
		this.sequence = sequence;
	}

	@Override
	public byte getBaseAt(String seqName, long position) {
		return sequence[(int) position];
	}

	@Override
	public void copyBases(String sequenceName, long from, int len, byte[] dest)
			throws IOException {
		System.arraycopy(sequence, (int) from, dest, 0, len) ;
	}

}
