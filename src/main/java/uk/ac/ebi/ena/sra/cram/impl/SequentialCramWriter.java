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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.zip.GZIPOutputStream;

import javax.swing.tree.DefaultMutableTreeNode;

import uk.ac.ebi.ena.sra.cram.CramException;
import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.encoding.BitCodec;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.compression.CramCompressionException;
import uk.ac.ebi.ena.sra.cram.impl.RecordCodecFactory.CodecStats;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.ExposedByteArrayOutputStream;

public class SequentialCramWriter {
	private RecordCodecFactory recordCodecFactory = new RecordCodecFactory();
	private CramRecordBlockWriter blockHeaderWriter;
	private BitOutputStream bos;
	private BitCodec<CramRecord> recordCodec;
	private OutputStream os;
	private SequenceBaseProvider referenceBaseProvider;

	private ExposedByteArrayOutputStream checkBAOS = new ExposedByteArrayOutputStream(1024);
	private BitOutputStream checkBOS = new DefaultBitOutputStream(checkBAOS);
	private ByteArrayInputStream checkBAIS = new ByteArrayInputStream(checkBAOS.getBuffer());
	private DefaultBitInputStream checkBIS = new DefaultBitInputStream(checkBAIS);
	private BitCodec<CramRecord> checkRecordReadCodec;
	private BitCodec<CramRecord> checkRecordWriteCodec;
	private RestoreBases restoreBases;

	private DefaultMutableTreeNode rootNode;
	private final CramHeader header;

	public int gzippedBlockHeaderBytes = 0;

	public SequentialCramWriter(OutputStream os, SequenceBaseProvider referenceBaseProvider, CramHeader header) {
		this.os = os;
		this.referenceBaseProvider = referenceBaseProvider;
		this.header = header;
		bos = new DefaultBitOutputStream(os);
		// bos = new DebuggingBitOuputStream(System.out, '\n', new
		// DefaultBitOutputStream(os));
		blockHeaderWriter = new CramRecordBlockWriter(os);
	}

	public long write(CramRecordBlock block) throws IOException, CramCompressionException {
		bos.flush();

		restoreBases = new RestoreBases(referenceBaseProvider, block.getSequenceName());

		rootNode = recordCodecFactory.buildCodecTree(header, block, referenceBaseProvider);
		checkRecordReadCodec = (BitCodec<CramRecord>) rootNode.getUserObject();

		rootNode = recordCodecFactory.buildCodecTree(header, block, referenceBaseProvider);
		checkRecordWriteCodec = (BitCodec<CramRecord>) rootNode.getUserObject();

		rootNode = recordCodecFactory.buildCodecTree(header, block, referenceBaseProvider);
		recordCodec = (BitCodec<CramRecord>) rootNode.getUserObject();

		// recordCodec = recordCodecFactory.createRecordCodec(header, block,
		// referenceBaseProvider) ;

		long bits = 0;
		bits = blockHeaderWriter.write(block);

		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		GZIPOutputStream gos = new GZIPOutputStream(baos);
		CramRecordBlockWriter bw = new CramRecordBlockWriter(gos);
		bw.write(block);
		gos.close();
		gzippedBlockHeaderBytes = baos.size();

		return bits;
	}

	public long write(CramRecord record) throws IOException {
		return recordCodec.write(bos, record);
	}

	public long writeAndCheck(CramRecord record) throws IOException, CramException {
		checkBAOS.reset();
		checkRecordWriteCodec.write(checkBOS, record);
		checkBOS.flush();

		checkBAIS.reset();
		checkBIS.reset();
		CramRecord checkRecord = checkRecordReadCodec.read(checkBIS);
		// the codec does not handle alignment start of the first sequence in
		// the block:
		checkRecord.setAlignmentStart(record.getAlignmentStart());
		if (checkRecord.isReadMapped())
			restoreBases.restoreReadBases(checkRecord);

		if (!checkRecord.equals(record)) {
			System.err.println(record.toString());
			System.err.println(checkRecord.toString());
			checkBAIS.reset();
			checkRecord = checkRecordReadCodec.read(checkBIS);

			throw new CramException("Round trip check failed.");
		}

		return write(record);
	}

	public void flush() throws IOException {
		bos.flush();
		os.flush();
	}

	public void dump() {
		RecordCodecFactory.dump(rootNode);
	}

	public CodecStats getCodecStats() {
		return RecordCodecFactory.getCodecStats(rootNode);
	}
}
