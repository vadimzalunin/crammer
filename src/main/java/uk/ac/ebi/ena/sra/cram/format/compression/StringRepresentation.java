package uk.ac.ebi.ena.sra.cram.format.compression;

public class StringRepresentation {

	public static String[] parse(String spec) {
		return spec.split(",");
	}

	public static String toString(Object... values) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < values.length; i++) {
			Object value = values[i];

			if (value instanceof Boolean)
				sb.append((Boolean) value ? "1" : "0");
			else
				sb.append(value);

			if (i < values.length - 1)
				sb.append(",");
		}
		return sb.toString();
	}

	public static boolean toBoolean(String value)
			throws CramCompressionException {
		if ("1".equals(value))
			return true;
		if ("0".equals(value))
			return false;
		throw new CramCompressionException("Boolean (1 or 0) value expected: "
				+ value);
	}

	public static long toLong(String value) throws CramCompressionException {
		try {
			return Long.valueOf(value);
		} catch (NumberFormatException e) {
			throw new CramCompressionException(
					"Expecting long value: " + value, e);
		}
	}

	public static int toInt(String value) throws CramCompressionException {
		try {
			return Integer.valueOf(value);
		} catch (NumberFormatException e) {
			throw new CramCompressionException(
					"Expecting long value: " + value, e);
		}
	}
	

	public static String numberEncodingToString(NumberEncoding encoding)
			throws CramCompressionException {
		switch (encoding) {
		case UNARY:
			return "UN";
		case BETA:
			return "BE";
		case GAMMA:
			return "GA";
		case GOLOMB:
			return "GO";
		case GOLOMB_RICE:
			return "GR";
		case SUBEXP:
			return "SU";

		default:
			throw new CramCompressionException(
					"Unsupported number compression: " + encoding);
		}
	}

	public static NumberEncoding stringToNumberEncoding(String string)
			throws CramCompressionException {
		if (string.length() != 2)
			throw new CramCompressionException(
					"Expecting a two-letter string: " + string);

		for (NumberEncoding encoding : NumberEncoding.values())
			if (string.equals(numberEncodingToString(encoding)))
				return encoding;

		throw new CramCompressionException("Unkown number encoding: " + string);
	}
}
