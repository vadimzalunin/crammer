package uk.ac.ebi.ena.sra.cram.format.compression;

import uk.ac.ebi.ena.sra.cram.encoding.GolombCodec;

class GolombCodecStub extends GolombCodec implements NumberCodecStub {

	public GolombCodecStub() {
		super (2) ;
	}

	@Override
	public NumberEncoding getEncoding() {
		return NumberEncoding.GOLOMB;
	}

	@Override
	public String getStringRepresentation() {
		return String.format("%d,%d,%d", getM(), getOffset(),
				isQuotientBit() ? 1 : 0);
	}

	@Override
	public void initFromString(String spec) throws CramCompressionException {
		String[] params = StringRepresentation.parse(spec);
		switch (params.length) {
		case 1:
			int m = StringRepresentation.toInt(params[0]);
			setM(m);
			break ;
		case 3:
			m = StringRepresentation.toInt(params[0]);
			setM(m);
			long offset = StringRepresentation.toLong(params[1]);
			setOffset(offset);
			boolean quotientBit = StringRepresentation.toBoolean(params[2]);
			setQuotientBit(quotientBit);
			break ;
		default:
			throw new CramCompressionException(
					"Not supported number of parameters to golomb codec: "
							+ params.length);
		}
	}

}
