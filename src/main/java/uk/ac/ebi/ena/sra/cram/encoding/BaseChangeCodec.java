package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class BaseChangeCodec implements BitCodec<BaseChange> {

	@Override
	public BaseChange read(BitInputStream bis) throws IOException {
		return new BaseChange(bis.readBits(2));
	}

	@Override
	public long write(BitOutputStream bis, BaseChange baseChange) throws IOException {
		bis.write(baseChange.getChange(), 2);
		return 2;
	}

	@Override
	public long numberOfBits(BaseChange baseChange) {
		return 2;
	}

}
