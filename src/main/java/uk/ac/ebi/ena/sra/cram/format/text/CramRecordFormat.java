package uk.ac.ebi.ena.sra.cram.format.text;

import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.util.List;

import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;

public class CramRecordFormat {
	public static final byte FIELD_SEPARATOR = '\t';
	public static final String STRING_FIELD_SEPARATOR = new String(
			new byte[] { FIELD_SEPARATOR });
	public static final byte NO_VALUE = '*';
	private ByteBuffer buf = ByteBuffer.allocate(1024 * 1024);
	private boolean eolFound = false;
	private boolean eofFound = false;
	private ReadFeaturesFormat rfFormat = new DefaultReadFeaturesFormat();

	public static void main(String[] args) throws IOException {

		CramRecordFormat format = new CramRecordFormat();

		CramRecord record;
		record = format.fromString("123	*	POS	*	*	*");
		System.out.println(record.toString());
		System.out.println(format.writeRecord(record));

		record = format
				.fromString("124	36	POS	M2wwxyzT$M2IAC.M5D2A!C!M13A!C!	*	*");
		System.out.println(record.toString());
		System.out.println(format.writeRecord(record));

		record = format
				.fromString("125	36	POS	M2wwxyzT$M2IAC.M5D2A!C!M13A!C!	GTGCGGATGCTCTCCTCCAGTTTGGGCTCGTGGTGTGTGTCCAGCAGGGACTGG	BBBBBBB=BBBBBBBBBBBBBBBB?BBBBB?BB?BBBBB?BBBBB?BBB??BBB");
		System.out.println(record.toString());
		System.out.println(format.writeRecord(record));
	}

	public String writeRecord(CramRecord record) throws IOException {
		StringBuilder sb = new StringBuilder();
		if (record.getAlignmentStart() > 0)
			sb.append(record.getAlignmentStart());
		else
			sb.appendCodePoint(NO_VALUE);

		sb.appendCodePoint(FIELD_SEPARATOR);

		if (record.getReadLength() > 0)
			sb.append(record.getReadLength());
		else
			sb.appendCodePoint(NO_VALUE);

		sb.appendCodePoint(FIELD_SEPARATOR);

		if (record.getRecordsToNextFragment() > 0)
			sb.append(record.getRecordsToNextFragment());
		else
			sb.appendCodePoint(NO_VALUE);

		sb.appendCodePoint(FIELD_SEPARATOR);

		if (record.isNegativeStrand())
			sb.append("NEG");
		else
			sb.append("POS");

		sb.appendCodePoint(FIELD_SEPARATOR);

		if (record.getReadFeatures() == null
				|| record.getReadFeatures().isEmpty())
			sb.appendCodePoint(NO_VALUE);
		else
			rfFormat.addFeaturesToStringBuilder(record.getReadFeatures(),
					record.getReadLength() > 0 ? (int) record.getReadLength()
							: 0, sb);

		sb.appendCodePoint(FIELD_SEPARATOR);

		if (record.getReadBases() == null)
			sb.appendCodePoint(NO_VALUE);
		else
			sb.append(new String(record.getReadBases()));

		sb.appendCodePoint(FIELD_SEPARATOR);

		if (record.getQualityScores() == null)
			sb.appendCodePoint(NO_VALUE);
		else
			sb.append(new String(record.getQualityScores()));

		return sb.toString();
	}

	public CramRecord fromString(String string) {
		CramRecord record = new CramRecord();

		String[] chunks = string.split(STRING_FIELD_SEPARATOR);

		int i = 0;
		String chunk;

		chunk = chunks[i++];
		if (isNotEmpty(chunk)) {
			record.setAlignmentStart(toLong(chunk));
			record.setReadMapped(true);
		} else
			record.setReadMapped(false);

		chunk = chunks[i++];
		if (isNotEmpty(chunk))
			record.setReadLength(toLong(chunk));

		chunk = chunks[i++];
		if (isNotEmpty(chunk)) {
			record.setRecordsToNextFragment(toLong(chunk));
			record.setLastFragment(false);
		} else
			record.setLastFragment(true);

		chunk = chunks[i++];
		if (chunk.equals("NEG"))
			record.setNegativeStrand(true);
		else
			record.setNegativeStrand(false);

		chunk = chunks[i++];
		if (isNotEmpty(chunk)) {
			List<ReadFeature> features = rfFormat.asReadFeatureList(new String(
					chunk));
			record.setReadFeatures(features);
			record.setPerfectMatch(false);
		} else
			record.setPerfectMatch(record.isReadMapped());

		chunk = chunks[i++];
		if (isNotEmpty(chunk))
			record.setReadBases(chunk.getBytes());

		chunk = chunks[i++];
		if (isNotEmpty(chunk))
			record.setQualityScores(chunk.getBytes());

		return record;
	}

	private CramRecord readRecord(InputStream is) throws IOException {
		CramRecord record = new CramRecord();
		eolFound = false;
		eofFound = false;

		byte[] start = readLineUntil(is, FIELD_SEPARATOR);
		if (isNotEmpty(start))
			record.setAlignmentStart(toLong(start));

		byte[] len = readLineUntil(is, FIELD_SEPARATOR);
		if (isNotEmpty(len))
			record.setReadLength(toLong(len));

		byte[] neg = readLineUntil(is, FIELD_SEPARATOR);
		if (new String(neg).equals("NEG"))
			record.setNegativeStrand(true);
		else
			record.setNegativeStrand(false);

		byte[] superCigar = readLineUntil(is, FIELD_SEPARATOR);
		if (isNotEmpty(superCigar)) {
			List<ReadFeature> features = rfFormat.asReadFeatureList(new String(
					superCigar));
			record.setReadFeatures(features);
		}

		byte[] bases = readLineUntil(is, FIELD_SEPARATOR);
		if (isNotEmpty(bases))
			record.setReadBases(bases);

		byte[] scores = readLineUntil(is, FIELD_SEPARATOR);
		if (isNotEmpty(scores))
			record.setQualityScores(scores);

		return record;
	}

	private static final boolean isNotEmpty(String s) {
		return s.length() != 1 || NO_VALUE != s.getBytes()[0];
	}

	private static final boolean isNotEmpty(byte[] bytes) {
		return bytes.length != 1 || bytes[0] != '*';
	}

	private static final long toLong(String s) {
		return Long.valueOf(s);
	}

	private static final long toLong(byte[] bytes) {
		return Long.valueOf(new String(bytes));
	}

	private byte[] readLineUntil(InputStream is, byte stopByte)
			throws IOException {
		if (eolFound)
			throw new RuntimeException("Premature end of line.");
		if (eofFound)
			throw new RuntimeException("Premature end of file.");

		byte b = 0;
		buf.clear();
		while ((b = (byte) is.read()) != -1 && b != stopByte && b != '\n')
			buf.put(b);

		eolFound = b == '\n';
		eofFound = b == -1;

		buf.flip();
		byte[] bytes = new byte[buf.limit()];
		buf.get(bytes);
		return bytes;
	}
}
