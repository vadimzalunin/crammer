package uk.ac.ebi.ena.sra.cram.format.compression;

import uk.ac.ebi.ena.sra.cram.encoding.BetaCodec;

class NullCodecStub extends BetaCodec implements NumberCodecStub {

	@Override
	public EncodingAlgorithm getEncoding() {
		return EncodingAlgorithm.NULL;
	}

	@Override
	public String getStringRepresentation() {
		return "";
	}

	@Override
	public void initFromString(String spec) throws CramCompressionException {
		if (spec != null && spec.length() > 1)
			throw new CramCompressionException(
					"Null codec does not expect any parameteres: " + spec);
	}

	@Override
	public Object[] getParameters() {
		return new Object[] {};
	}

	@Override
	public void setParameters(Object[] params) {
	}
}
