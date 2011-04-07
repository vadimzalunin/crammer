package uk.ac.ebi.ena.sra.cram.format.compression;

import uk.ac.ebi.ena.sra.cram.encoding.GammaCodec;

class GammaCodecStub extends GammaCodec implements NumberCodecStub {

	public GammaCodecStub() {
		super();
	}

	@Override
	public NumberEncoding getEncoding() {
		return NumberEncoding.GAMMA;
	}

	@Override
	public String getStringRepresentation() {
		return String.format("%d,%d", getOffset(),
				isLenCodingBit() ? 1 : 0);
	}

	@Override
	public void initFromString(String spec) throws CramCompressionException {
		String[] params = StringRepresentation.parse(spec);
		switch (params.length) {
		case 0:
			break;
		case 2:
			boolean lenCodingBit = StringRepresentation.toBoolean(params[0]);
			setLenCodingBit(lenCodingBit) ;
			long offset = StringRepresentation.toLong(params[1]);
			setOffset(offset);
			break ;
		default:
			throw new CramCompressionException(
					"Not supported number of parameters to unary codec: "
							+ params.length);
		}
	}

}
