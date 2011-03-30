package uk.ac.ebi.ena.sra.cram.impl;

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

}
