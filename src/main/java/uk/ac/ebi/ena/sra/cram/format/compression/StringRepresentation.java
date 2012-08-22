/*******************************************************************************
 * Copyright 2012 EMBL-EBI, Hinxton outstation
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *   http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/
package uk.ac.ebi.ena.sra.cram.format.compression;

class StringRepresentation {

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

	public static boolean toBoolean(String value) throws CramCompressionException {
		if ("1".equals(value))
			return true;
		if ("0".equals(value))
			return false;
		throw new CramCompressionException("Boolean (1 or 0) value expected: " + value);
	}

	public static long toLong(String value) throws CramCompressionException {
		try {
			return Long.valueOf(value);
		} catch (NumberFormatException e) {
			throw new CramCompressionException("Expecting long value: " + value, e);
		}
	}

	public static int toInt(String value) throws CramCompressionException {
		try {
			return Integer.valueOf(value);
		} catch (NumberFormatException e) {
			throw new CramCompressionException("Expecting long value: " + value, e);
		}
	}

	public static String encodingToString(EncodingAlgorithm encoding) throws CramCompressionException {
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
			throw new CramCompressionException("Unsupported number compression: " + encoding);
		}
	}

	public static EncodingAlgorithm stringToNumberEncoding(String string) throws CramCompressionException {
		if (string.length() != 2)
			throw new CramCompressionException("Expecting a two-letter string: " + string);

		for (EncodingAlgorithm encoding : EncodingAlgorithm.values())
			if (string.equals(encodingToString(encoding)))
				return encoding;

		throw new CramCompressionException("Unkown number encoding: " + string);
	}
}
