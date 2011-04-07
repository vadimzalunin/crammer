package uk.ac.ebi.ena.sra.cram.format.compression;

import uk.ac.ebi.ena.sra.cram.encoding.BetaCodec;

class BetaCodecStub extends BetaCodec implements NumberCodecStub {

	@Override
	public NumberEncoding getEncoding() {
		return NumberEncoding.BETA;
	}

	@Override
	public String getStringRepresentation() {
		return String.format("%d,%d", getOffset(), getReadNofBits());
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
			int readNofBits = StringRepresentation.toInt(params[1]);
			setReadNofBits(readNofBits);
			break ;
		default:
			throw new CramCompressionException(
					"Not supported number of parameters to beta codec: "
							+ params.length);
		}
	}

}
