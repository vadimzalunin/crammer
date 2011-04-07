package uk.ac.ebi.ena.sra.cram.format.compression;

import uk.ac.ebi.ena.sra.cram.encoding.UnaryCodec;

class UnaryCodecStub extends UnaryCodec implements NumberCodecStub {

	public UnaryCodecStub() {
		super();
	}

	@Override
	public NumberEncoding getEncoding() {
		return NumberEncoding.UNARY;
	}

	@Override
	public String getStringRepresentation() {
		return String.format("%d,%d", getOffset(), isStopBit() ? 1 : 0);
	}

	@Override
	public void initFromString(String spec) throws CramCompressionException {
		String[] params = StringRepresentation.parse(spec);
		switch (params.length) {
		case 0:
			break;
		case 2:
			boolean stopBit = StringRepresentation.toBoolean(params[0]);
			setStopBit(stopBit);
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
