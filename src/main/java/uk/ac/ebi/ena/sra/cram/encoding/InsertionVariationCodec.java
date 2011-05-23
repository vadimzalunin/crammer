package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.NullBitOutputStream;

public class InsertionVariationCodec implements BitCodec<InsertionVariation> {
	public BitCodec<byte[]> insertBasesCodec;

	@Override
	public InsertionVariation read(BitInputStream bis) throws IOException {
		// position is not read here because we need to keep track of previous
		// values read from the codec. See ReadFeatureCodec.
		long position = -1L;
		byte[] insertion = insertBasesCodec.read(bis);

		InsertionVariation v = new InsertionVariation();
		v.setPosition((int) position);
		v.setSequence(insertion);
		return v;
	}

	@Override
	public long write(BitOutputStream bos, InsertionVariation v)
			throws IOException {
		long len = 0L;

		long seqLen = 0L ;
		seqLen -=len ;
		len += insertBasesCodec.write(bos, v.getSequence());
		seqLen += len ;
//		System.out.println(new String(v.getSequence()) + " encoded in " +  seqLen + " bits.");

		return len;
	}

	@Override
	public long numberOfBits(InsertionVariation v) {
		try {
			return write(NullBitOutputStream.INSTANCE, v);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

}
