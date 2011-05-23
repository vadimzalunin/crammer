package uk.ac.ebi.ena.sra.cram.format.compression;

import uk.ac.ebi.ena.sra.cram.encoding.GammaCodec;

class GammaCodecStub extends GammaCodec implements NumberCodecStub {

	public GammaCodecStub() {
		super();
	}

	@Override
	public EncodingAlgorithm getEncoding() {
		return EncodingAlgorithm.GAMMA;
	}

	@Override
	public String getStringRepresentation() {
		return String.format("%d,%d", getOffset(), isLenCodingBit() ? 1 : 0);
	}

	@Override
	public void initFromString(String spec) throws CramCompressionException {
		String[] params = StringRepresentation.parse(spec);
		switch (params.length) {
		case 0:
			break;
		case 2:
			long offset = StringRepresentation.toLong(params[0]);
			setOffset(offset);
			boolean lenCodingBit = StringRepresentation.toBoolean(params[1]);
			setLenCodingBit(lenCodingBit);
			break;
		default:
			throw new CramCompressionException(
					"Not supported number of parameters to unary codec: "
							+ params.length);
		}
	}

	@Override
	public Object[] getParameters() {
		return new Object[] { getOffset(), isLenCodingBit() };
	}

	@Override
	public void setParameters(Object[] params) {
		setOffset((Long) params[0]);
		setLenCodingBit((Boolean) params[1]);
	}
}
