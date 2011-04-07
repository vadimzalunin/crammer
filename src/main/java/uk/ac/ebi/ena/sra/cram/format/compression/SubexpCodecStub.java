package uk.ac.ebi.ena.sra.cram.format.compression;

import uk.ac.ebi.ena.sra.cram.encoding.SubexpCodec;

class SubexpCodecStub extends SubexpCodec implements NumberCodecStub {

	public SubexpCodecStub() {
		super(1);
	}

	@Override
	public NumberEncoding getEncoding() {
		return NumberEncoding.SUBEXP;
	}

	@Override
	public String getStringRepresentation() {
		return String.format("%d,%d,%d", getK(), getOffset(), isUnaryBit() ? 1
				: 0);
	}

	@Override
	public void initFromString(String spec) throws CramCompressionException {
		String[] params = StringRepresentation.parse(spec);
		switch (params.length) {
		case 1:
			int k = StringRepresentation.toInt(params[0]);
			setK(k);
			break;
		case 3:
			k = StringRepresentation.toInt(params[0]);
			setK(k);
			long offset = StringRepresentation.toLong(params[1]);
			setOffset(offset);
			boolean quotientBit = StringRepresentation.toBoolean(params[2]);
			setUnaryBit(quotientBit);
			break ;
		default:
			throw new CramCompressionException(
					"Not supported number of parameters to golomb-rice codec: "
							+ params.length);
		}
	}

}
