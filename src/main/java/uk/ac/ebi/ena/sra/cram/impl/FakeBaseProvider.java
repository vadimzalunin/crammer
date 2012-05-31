package uk.ac.ebi.ena.sra.cram.impl;

import java.io.IOException;
import java.util.Arrays;

import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;

public class FakeBaseProvider implements SequenceBaseProvider {
	private byte base = 'N' ; 

	@Override
	public byte getBaseAt(String seqName, long position) {
		return base ;
	}

	@Override
	public void copyBases(String sequenceName, long from, int len, byte[] dest) throws IOException {
		Arrays.fill(dest, 0, len, base) ;
	}

}
