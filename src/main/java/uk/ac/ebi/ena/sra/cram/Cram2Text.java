package uk.ac.ebi.ena.sra.cram;

import java.io.BufferedInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Scanner;

import uk.ac.ebi.ena.sra.cram.encoding.JavaBinaryCramReader;
import uk.ac.ebi.ena.sra.cram.format.CramFormatException;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.impl.ByteArraySequenceBaseProvider;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

public class Cram2Text {

	public static void main(String[] args) throws IOException,
			CramFormatException {
		Params params = new Params();
		JCommander jc = new JCommander(params);
		jc.parse(args);

		Scanner scanner = new Scanner(new FileInputStream(params.reference));
		ByteArrayOutputStream refSeqBAOS = new ByteArrayOutputStream();
		while (scanner.hasNextLine()) {
			String line = scanner.nextLine();
			if (line.startsWith(">"))
				continue;
			refSeqBAOS.write(line.getBytes());
		}

		ByteArraySequenceBaseProvider provider = new ByteArraySequenceBaseProvider(
				refSeqBAOS.toByteArray());
		refSeqBAOS.close();
		refSeqBAOS = null;

		InputStream is = new BufferedInputStream(new FileInputStream(
				params.cramFile));
		JavaBinaryCramReader reader = new JavaBinaryCramReader(provider,
				new DataInputStream(is));

		long seqRecordCounter = 1;
		while (reader.hasNext()) {
			CramRecord cramRecord = reader.next();
			System.out.println(Utils.toString(String.valueOf(seqRecordCounter),
					Utils.restoreBases(cramRecord, provider, reader
							.getCurrentBlock().getSequenceName()), cramRecord));
			if (params.maxRecords > -1
					&& seqRecordCounter++ > params.maxRecords)
				break;
		}
	}

	static class Params {
		@Parameter(names = { "--cram-file" })
		File cramFile;

		@Parameter(names = { "--max-records" })
		int maxRecords = -1;

		@Parameter(names = { "--reference" }, required = true)
		String reference;
	}
}
