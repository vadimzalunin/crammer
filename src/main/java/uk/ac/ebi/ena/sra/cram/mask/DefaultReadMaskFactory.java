package uk.ac.ebi.ena.sra.cram.mask;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.regex.Pattern;

public class DefaultReadMaskFactory implements ReadMaskFactory<String> {

	private static ReadMaskFactory<String> detectReadMaskFormat(InputStream is)
			throws IOException {
		BufferedReader bis = new BufferedReader(new InputStreamReader(is));
		String line = null;

		Pattern intMaskPattern = Pattern.compile("^[0-9\\s]+$");
		Pattern fastaMaskPattern = Pattern.compile("^["
				+ FastaByteArrayMaskFactory.DEFAULT_MASK_BYTE
				+ FastaByteArrayMaskFactory.DEFAULT_NON_MASK_BYTE + "]+$");
		while ((line = bis.readLine()) != null) {
			if (line.length() == 0)
				continue;
			boolean intFormatMatches = intMaskPattern.matcher(line).matches();
			boolean fastaFormatMatches = fastaMaskPattern.matcher(line)
					.matches();

			if (intFormatMatches && fastaFormatMatches)
				continue;

			if (intFormatMatches)
				return new IntegerListMaskFactory();

			return new FastaByteArrayMaskFactory();
		}
		return null;
	}

	@Override
	public PositionMask createMask(String line) throws ReadMaskFormatException {
		return null;
	}

}
