package uk.ac.ebi.ena.sra.cram;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;
import uk.ac.ebi.ena.sra.cram.format.InsertBase;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.compression.CramCompressionException;
import uk.ac.ebi.ena.sra.cram.impl.ByteArraySequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.impl.CramHeaderIO;
import uk.ac.ebi.ena.sra.cram.impl.SequentialCramWriter;
import uk.ac.ebi.ena.sra.cram.stats.CramStats;

public class CramExample {

	public static void main(String[] args) throws IOException,
			CramCompressionException {
		// assume there is only one reference sequence:
		String refSeqName = "1";
		byte[] referenceSequenceBases = "AAAAAAAAAAAAAAAAAAAA".getBytes();

		OutputStream os = new BufferedOutputStream(new FileOutputStream(
				new File(args[0])));
		
		CramHeader header = new CramHeader();
		// prepare header here:
		header.setReferenceSequences(Arrays.asList(new CramReferenceSequence(
				refSeqName, referenceSequenceBases.length)));
		// ...
		CramHeaderIO.write(header, os);
		
		SequenceBaseProvider provider = new ByteArraySequenceBaseProvider(
				referenceSequenceBases);
		SequentialCramWriter writer = new SequentialCramWriter(os,
				provider, header);

		// assume we have 10 blocks to write:
		int nofBlocks = 10;
		for (int blockIndex = 0; blockIndex < nofBlocks; blockIndex++) {
			CramRecordBlock block = new CramRecordBlock();
			// configure block params here: 
			block.setSequenceName(refSeqName) ;
			// ...
			
			CramStats stats = new CramStats(header, null);

			// assume there are 100k records in the block:
			int recordsInBlock = 100000;
			List<CramRecord> blockRecords = new ArrayList<CramRecord>(recordsInBlock) ;
			for (int recordIndex = 0; recordIndex < recordsInBlock; recordIndex++) {
				CramRecord record = new CramRecord();
				// set record fields here:
				record.setReadFeatures(Arrays
						.asList(new ReadFeature[] { new InsertBase(1,
								(byte) 'A') }));
				// ...

				stats.addRecord(record);
				blockRecords.add(record);
			}

			stats.adjustBlock(block);
			writer.write(block);
			
			for (CramRecord record : block.getRecords())
				writer.write(record);

		}
		writer.flush();
	}
}
