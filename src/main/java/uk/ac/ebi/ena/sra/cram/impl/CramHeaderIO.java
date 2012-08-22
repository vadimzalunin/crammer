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
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMTag;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.format.CramFormatException;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramHeaderRecord;
import uk.ac.ebi.ena.sra.cram.format.CramReadGroup;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;
import uk.ac.ebi.ena.sra.cram.format.ReadAnnotation;

public class CramHeaderIO {
	private static final String MAGIC = "CRAM";

	public static void write(CramHeader header, OutputStream os) throws IOException {
		DataOutputStream headerDOS = new DataOutputStream(os);
		headerDOS.writeUTF(MAGIC);
		headerDOS.writeUTF(header.getVersion());

		headerDOS.writeInt(header.getReferenceSequences().size());
		for (CramReferenceSequence s : header.getReferenceSequences()) {
			headerDOS.writeUTF(s.getName());
			headerDOS.writeInt(s.getLength());
		}

		if (header.getReadAnnotations() == null || header.getReadAnnotations().isEmpty())
			headerDOS.writeInt(0);
		else {
			headerDOS.writeInt(header.getReadAnnotations().size());
			for (ReadAnnotation ra : header.getReadAnnotations())
				headerDOS.writeUTF(ra.getKey());
		}

		if (header.getReadGroups() == null || header.getReadGroups().isEmpty())
			headerDOS.writeInt(0);
		else {
			headerDOS.writeInt(header.getReadGroups().size());
			for (CramReadGroup rg : header.getReadGroups()) {
				if (rg.getId() == null)
					headerDOS.writeUTF("");
				else
					headerDOS.writeUTF(rg.getId());

				if (rg.getSample() == null)
					headerDOS.writeUTF("");
				else
					headerDOS.writeUTF(rg.getSample());
			}
		}

		List<CramHeaderRecord> records = header.getRecords();
		if (records == null || records.isEmpty()) {
			headerDOS.writeInt(0);
		} else {
			headerDOS.writeInt(records.size());
			for (CramHeaderRecord record : records) {
				byte[] tagBytes = record.getTag().getBytes();
				if (tagBytes.length != 2)
					throw new RuntimeException("Header tag must be 2 bytes long: " + record.getTag());
				headerDOS.write(tagBytes);

				Set<String> keySet = record.getKeySet();
				headerDOS.writeInt(keySet.size());
				for (String key : keySet) {
					byte[] keyBytes = key.getBytes();
					if (keyBytes.length != 2)
						throw new RuntimeException("Header tag key must be 2 bytes long: " + key);
					headerDOS.write(keyBytes);

					String value = record.getValue(key);
					headerDOS.writeUTF(value);
				}
			}
		}

		headerDOS.flush();
	}

	public static CramHeader read(InputStream is) throws IOException, CramFormatException {
		CramHeader header = new CramHeader();
		DataInputStream dis = null;
		if (is instanceof DataInputStream)
			dis = (DataInputStream) is;
		else
			dis = new DataInputStream(is);

		String magick = dis.readUTF();
		if (!MAGIC.equals(magick))
			throw new RuntimeException("Not recognized as CRAM format.");

		String version = dis.readUTF();
		if (!Utils.isCompatible(version))
			throw new CramFormatException("incompatble CRAM data version " + version);

		header.setVersion(version);
		int seqCount = dis.readInt();
		ArrayList<CramReferenceSequence> seqList = new ArrayList<CramReferenceSequence>();
		for (int i = 0; i < seqCount; i++) {
			CramReferenceSequence cramSeq = new CramReferenceSequence();
			cramSeq.setName(dis.readUTF());
			cramSeq.setLength(dis.readInt());
			seqList.add(cramSeq);
		}
		header.setReferenceSequences(seqList);

		int annSize = dis.readInt();
		List<ReadAnnotation> annotations = new ArrayList<ReadAnnotation>(annSize);
		for (int i = 0; i < annSize; i++) {
			String annKey = dis.readUTF();
			annotations.add(new ReadAnnotation(annKey));
		}
		header.setReadAnnotations(annotations);

		int rgSize = dis.readInt();
		List<CramReadGroup> rgList = new ArrayList<CramReadGroup>(rgSize);
		for (int i = 0; i < rgSize; i++) {
			String id = dis.readUTF();
			if (id.length() == 0)
				id = null;
			String sample = dis.readUTF();
			if (sample.length() == 0)
				sample = null;
			rgList.add(new CramReadGroup(id, sample));
		}
		header.setReadGroups(rgList);

		int recordCount = dis.readInt();
		List<CramHeaderRecord> records = header.getRecords();
		Set<String> seqRecordSet = new TreeSet<String>();
		Set<String> rgRecordSet = new TreeSet<String>();
		for (int recordIndex = 0; recordIndex < recordCount; recordIndex++) {
			byte[] tagBytes = new byte[2];
			dis.readFully(tagBytes);
			String tag = new String(tagBytes);
			CramHeaderRecord record = new CramHeaderRecord(tag);
			records.add(record);

			int recordSize = dis.readInt();
			for (int keyIndex = 0; keyIndex < recordSize; keyIndex++) {
				byte[] keyBytes = new byte[2];
				dis.readFully(keyBytes);
				String key = new String(keyBytes);
				String value = dis.readUTF();
				record.setValue(key, value);

				if (SAMTag.SQ.name().equals(tag) && SAMSequenceRecord.SEQUENCE_NAME_TAG.equals(key))
					seqRecordSet.add(value);

				if (SAMTag.RG.name().equals(tag) && SAMReadGroupRecord.READ_GROUP_ID_TAG.equals(key))
					seqRecordSet.add(value);

			}
		}

		// merge sequences, a workaround:
		if (header.getReferenceSequences() != null && !header.getReferenceSequences().isEmpty()) {

			for (CramReferenceSequence seq : header.getReferenceSequences()) {
				if (seqRecordSet.contains(seq.getName()))
					continue;
				CramHeaderRecord record = new CramHeaderRecord(SAMTag.SQ.name());
				record.setValue(SAMSequenceRecord.SEQUENCE_LENGTH_TAG, String.valueOf(seq.getLength()));
				record.setValue(SAMSequenceRecord.SEQUENCE_NAME_TAG, seq.getName());

				records.add(record);
			}

			for (CramReadGroup rg : header.getReadGroups()) {
				if (rg.getId() == null || seqRecordSet.contains(rg.getId()))
					continue;
				CramHeaderRecord record = new CramHeaderRecord(SAMTag.RG.name());
				record.setValue(SAMReadGroupRecord.READ_GROUP_ID_TAG, rg.getId());
				record.setValue(SAMReadGroupRecord.READ_GROUP_SAMPLE_TAG, rg.getSample());

				records.add(record);
			}
		}

		return header;
	}
}
