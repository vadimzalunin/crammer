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
package uk.ac.ebi.ena.sra.cram.impl;

import java.io.DataInputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SeekableFileStream;
import net.sf.samtools.util.SeekableStream;
import uk.ac.ebi.ena.sra.cram.CramException;
import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.encoding.BitCodec;
import uk.ac.ebi.ena.sra.cram.format.CramFormatException;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.compression.CramCompressionException;
import uk.ac.ebi.ena.sra.cram.index.CramIndex;
import uk.ac.ebi.ena.sra.cram.index.RecordPointer;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;

public class CramRandomAccessIterator implements CloseableIterator<CramRecord> {

	private SeekableStream seekableStream;
	private DataInputStream dis;
	private CramHeader header;
	private CramRecordBlock block;
	private boolean eof = false;

	private long recordCounter = 0;

	private RecordCodecFactory recordCodecFactory = new RecordCodecFactory();
	private SequenceBaseProvider referenceBaseProvider;
	private BitCodec<CramRecord> recordCodec;
	private DefaultBitInputStream bis;
	private RestoreBases restoreBases;
	private final RecordPointer start;
	private DataInputStream blockDIS;
	private ReferenceSequenceFile referenceSequenceFile;
	private ReferenceSequence referenceSequence;

	public static void main(String[] args) throws CramException, IOException {
		File cramFile = new File(
				"c:/temp/HG00096.mapped.illumina.mosaik.GBR.exome.20110411.chr20.cram");
		File cramIndexFile = new File(cramFile.getAbsolutePath() + ".crai");
		CramIndex index = CramIndex.fromFile(cramIndexFile);

		File refFile = new File("c:/temp/chr20.fa");
		ReferenceSequenceFile referenceSequenceFile = Utils
				.createIndexedFastaSequenceFile(refFile);

		byte[] refBases = Utils.getReferenceSequenceBases(referenceSequenceFile, "20");
		ByteArraySequenceBaseProvider provider = new ByteArraySequenceBaseProvider(
				refBases);

		Random random = new Random();
		long time1 = System.currentTimeMillis();
		for (int i = 0; i < 1000; i++) {
			long queryStart = random.nextInt(62905065);
			getRecordsStartingFrom(cramFile, index, referenceSequenceFile,
					"20", queryStart, 1000);
		}
		long time2 = System.currentTimeMillis();
		System.out.println(time2 - time1);
	}

	private static List<CramRecord> getRecordsStartingFrom(File cramFile,
			CramIndex index, ReferenceSequenceFile referenceSequenceFile,
			String seqName, long queryStart, int nofRecords)
			throws CramFormatException, CramCompressionException, IOException {

		SeekableFileStream stream = new SeekableFileStream(cramFile);

		List<CramRecord> records = new ArrayList<CramRecord>(nofRecords);

		RecordPointer pointer = index.findRecordPointerAt(seqName, queryStart);
		// System.out.println(pointer);
		CramRandomAccessIterator iterator = new CramRandomAccessIterator(
				stream, referenceSequenceFile, pointer);

		int counter = 0;
		while (iterator.hasNext()) {
			CramRecord record = iterator.next();
			if (record.getAlignmentStart() < queryStart)
				continue;
			records.add(record);
			if (++counter >= 10)
				break;
		}

		return records;
	}

	public CramRandomAccessIterator(SeekableStream stream,
			ReferenceSequenceFile referenceSequenceFile, RecordPointer start)
			throws IOException, CramFormatException, CramCompressionException {
		if (!Utils.isCRAM(stream)) throw new RuntimeException("Not a valid CRAM format.") ;
		
		this.seekableStream = stream;
		this.referenceSequenceFile = referenceSequenceFile;
		this.start = start;

		dis = new DataInputStream(stream);
		readHeader();
		seekableStream.seek(start.getBlockStart());
		readNextBlock();
		if (start.getBitOffset() == 0)
			blockDIS.skip(start.getByteOffset());
		else
			blockDIS.skip(start.getByteOffset() - 1);
		bis.readBits(start.getBitOffset());
		recordCounter = start.getRecordNumber();
	}

	private void readHeader() throws IOException, CramFormatException{
		header = CramHeaderIO.read(Utils.getNextChunk(dis));
	}

	private void readNextBlock() throws CramFormatException, IOException,
			CramCompressionException {
		blockDIS = Utils.getNextChunk(this.dis);
		if (blockDIS == null) {
			eof = true;
			return;
		}
		CramRecordBlockReader crbReader = new CramRecordBlockReader(blockDIS);
		block = crbReader.read();
		recordCounter = 0;

		if (block == null)
			eof = true;
		else {
			if (referenceSequence == null
					|| !block.getSequenceName().equals(
							referenceSequence.getName())) {
				if (referenceSequence == null)
					referenceSequence = referenceSequenceFile.getSequence(block
							.getSequenceName());

				referenceSequence = referenceSequenceFile.getSubsequenceAt(
						block.getSequenceName(), 1, referenceSequence.length());
				referenceBaseProvider = new ByteArraySequenceBaseProvider(
						referenceSequence.getBases());

			}

			recordCodec = recordCodecFactory.createRecordCodec(header, block,
					referenceBaseProvider);
			restoreBases = new RestoreBases(referenceBaseProvider,
					block.getSequenceName());
			bis = new DefaultBitInputStream(blockDIS);
		}
	}

	@Override
	public boolean hasNext() {
		if (eof)
			return false;

		if (block == null) {
			try {
				readNextBlock();
			} catch (Exception e) {
				throw new RuntimeException(e);
			}

			if (eof)
				return false;
		}

		return recordCounter < block.getRecordCount();
	}

	@Override
	public CramRecord next() {
		if (eof)
			throw new RuntimeException("End of stream.");

		if (block == null || recordCounter >= block.getRecordCount())
			throw new RuntimeException("Iterator depleted or invalid state.");

		CramRecord record;
		try {
			record = recordCodec.read(bis);
			record.setSequenceName(block.getSequenceName()) ;
			record.setAlignmentStart(record.getAlignmentStart()
					+ start.getAlignmentStart()
					- block.getFirstRecordPosition());
			recordCounter++;
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		if (record.isReadMapped())
			try {
				restoreBases.restoreReadBases(record);
			} catch (IOException e) {
				throw new RuntimeException(e);
			}

		return record;
	}

	@Override
	public void remove() {
		throw new RuntimeException("remove unsupported in CramRecordIterator.");
	}

	@Override
	public void close() {
		// TODO Auto-generated method stub
		
	}

}
