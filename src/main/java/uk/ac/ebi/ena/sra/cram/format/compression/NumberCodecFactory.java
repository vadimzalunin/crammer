package uk.ac.ebi.ena.sra.cram.format.compression;

public class NumberCodecFactory {

	public static NumberCodecStub createStub(EncodingAlgorithm encoding)
			throws CramCompressionException {
		switch (encoding) {
		case GOLOMB:
			return new GolombCodecStub();
		case UNARY:
			return new UnaryCodecStub();
		case GAMMA:
			return new GammaCodecStub();
		case GOLOMB_RICE:
			return new GolombRiceCodecStub();
		case SUBEXP:
			return new SubexpCodecStub();
		case BETA:
			return new BetaCodecStub();

		default:
			throw new CramCompressionException("Unsupported codec: " + encoding);
		}
	}
}