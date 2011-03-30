package uk.ac.ebi.ena.sra.cram;

import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.Scanner;

import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import uk.ac.ebi.ena.sra.cram.encoding.JavaBinaryCramReader;
import uk.ac.ebi.ena.sra.cram.format.CramFormatException;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.impl.ByteArraySequenceBaseProvider;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.FileConverter;

public class Cram2Bam {

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

		BAMFileWriter bamFileWriter = new BAMFileWriter(params.outputFile);
		SAMFileHeader header = new SAMFileHeader();
		bamFileWriter.setHeader(header);

		InputStream is = new FileInputStream(params.cramFile);
		DataInputStream dis = new DataInputStream(is);

		int recordCounter = 0;

		JavaBinaryCramReader reader = new JavaBinaryCramReader(provider, dis);
		while (reader.hasNext()) {
			CramRecord cramRecord = reader.next();
			byte[] refBases = Utils.restoreBases(cramRecord, provider, reader
					.getCurrentBlock().getSequenceName());

			SAMRecord samRecord = new SAMRecord(header);
			samRecord.setAlignmentStart((int) cramRecord.getAlignmentStart());
			samRecord.setReadBases(refBases);
			samRecord.setReadName(String.valueOf(recordCounter++));
			samRecord.setReadNegativeStrandFlag(cramRecord.isNegativeStrand());
			Arrays.fill(refBases, (byte) '!');
			samRecord.setBaseQualities(refBases);
			bamFileWriter.addAlignment(samRecord);

			if (recordCounter > params.maxRecordsPerSequence)
				break;
		}

		bamFileWriter.close();
	}

	@Parameters(commandDescription = "CRAM printing and conversion")
	static class Params {
		@Parameter(names = { "--input-cram" }, converter = FileConverter.class)
		File cramFile;

		@Parameter(names = { "--max-sequences" })
		int maxSequences = 0;

		@Parameter(names = { "--max-records-per-sequence" })
		long maxRecordsPerSequence = 0;

		@Parameter(names = { "--reference" }, required = true, converter = FileConverter.class)
		File reference;

		@Parameter(names = { "--output-file" }, converter = FileConverter.class)
		File outputFile;
	}
}
