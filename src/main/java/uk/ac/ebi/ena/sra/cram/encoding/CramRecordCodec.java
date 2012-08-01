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
package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.ReadAnnotation;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.ReadTag;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.NullBitOutputStream;

public class CramRecordCodec implements BitCodec<CramRecord> {
	public BitCodec<Long> inSeqPosCodec;
	public BitCodec<Long> recordsToNextFragmentCodec;
	public BitCodec<Long> readlengthCodec;
	public BitCodec<List<ReadFeature>> variationsCodec;
	public SequenceBaseProvider sequenceBaseProvider;

	public BitCodec<Byte> baseCodec;
	public ByteArrayBitCodec qualityCodec;

	public String sequenceName;
	public long prevPosInSeq = 1L;
	public long defaultReadLength = 0L;

	public BitCodec<ReadAnnotation> readAnnoCodec;
	public BitCodec<Integer> readGroupCodec;
	public BitCodec<Long> nextFragmentIDCodec;

	public BitCodec<Byte> mappingQualityCodec;

	public boolean storeMappedQualityScores = false;
	public BitCodec<Byte> heapByteCodec;

	public TreeMap<String, BitCodec<byte[]>> tagCodecMap;
	public BitCodec<String> tagKeyAndTypeCodec;
	public Map<String, BitCodec<Integer>> tagValueByteLenCodec;
	public BitCodec<Byte> tagCountCodec;

	public BitCodec<Byte> flagsCodec;

	private static Logger log = Logger.getLogger(CramRecordCodec.class);

	private static int debugRecordEndMarkerLen = 0;
	private static long debugRecordEndMarker = ~(-1 << (debugRecordEndMarkerLen / 2));

	public BitCodec<byte[]> readNameCodec;
	public boolean preserveReadNames = false;

	private List<Integer> tagLens = new ArrayList<Integer>();

	private List<String> deferredTags = new ArrayList<String>();

	private long[] timePoints = new long[10];
	private long time = 0;
	private int timeIndex = 0;
	

	private void timeTick() {
		long stop = System.nanoTime();
		timePoints[timeIndex++] += stop - time;
		time = stop;
	}
	
	long readTime = 0 ;

	@Override
	public CramRecord read(BitInputStream bis) throws IOException {
		long startMillis = System.currentTimeMillis() ;
		
		timeIndex = 0;
		time = System.nanoTime();

		recordCounter++;

		// long marker = bis.readLongBits(debugRecordEndMarkerLen);
		// if (marker != debugRecordEndMarker) {
		// throw new
		// RuntimeException("Debug marker for beginning of record not found.");
		// }

		CramRecord record = new CramRecord();

		byte b = flagsCodec.read(bis);
		record.setFlags(b);

		String readName = null;
		if (preserveReadNames)
			readName = new String(readNameCodec.read(bis));

		if (!record.isLastFragment()) {
			if (!record.detached) {
				record.setRecordsToNextFragment(recordsToNextFragmentCodec.read(bis));
			} else {
				CramRecord mate = new CramRecord();
				mate.setReadMapped(bis.readBit());
				mate.setNegativeStrand(bis.readBit());
				mate.setFirstInPair(bis.readBit());
				if (readName == null)
					readName = new String(readNameCodec.read(bis));
				mate.setReadName(readName);
				mate.setSequenceName(readZeroTerminatedString(heapByteCodec, bis));
				mate.setAlignmentStart(Long.valueOf(readZeroTerminatedString(heapByteCodec, bis)));
				record.insertSize = Integer.valueOf(readZeroTerminatedString(heapByteCodec, bis));

				mate.setFirstInPair(!record.isFirstInPair());
				if (record.isFirstInPair())
					record.next = mate;
				else
					record.previous = mate;
			}
		}
		record.setReadName(readName);

		timeTick();

		int readLen;
		if (bis.readBit())
			readLen = readlengthCodec.read(bis).intValue();
		else
			readLen = (int) defaultReadLength;
		record.setReadLength(readLen);

		record.setReadGroupID(readGroupCodec.read(bis));

		long position = prevPosInSeq + inSeqPosCodec.read(bis);
		prevPosInSeq = position;
		record.setAlignmentStart(position);

		timeTick();

		int tagCount = tagCountCodec.read(bis);
		deferredTags.clear();
		tagLens.clear();
		if (tagCount > 0)
			record.tags = new ArrayList<ReadTag>();
		for (int i = 0; i < tagCount; i++) {
			String tagKeyAndType = tagKeyAndTypeCodec.read(bis);

			BitCodec<byte[]> tagValueCodec = tagCodecMap.get(tagKeyAndType);
			if (tagValueCodec != null) {
				byte[] value = tagValueCodec.read(bis);

				ReadTag tag = ReadTag.deriveTypeFromKeyAndType(tagKeyAndType,
						ReadTag.restoreValueFromByteArray(tagKeyAndType.charAt(3), value));
				record.tags.add(tag);
			} else {
				BitCodec<Integer> lenCodec = tagValueByteLenCodec.get(tagKeyAndType);
				int tagValueByteLen = lenCodec.read(bis);
				tagLens.add(tagValueByteLen);
				deferredTags.add(tagKeyAndType);
			}
		}

		timeTick();

		boolean bisByteAligned = false;

		if (record.isReadMapped()) {
			boolean imperfectMatch = bis.readBit();
			if (imperfectMatch) {
				List<ReadFeature> features = variationsCodec.read(bis);
				record.setReadFeatures(features);
			}

			record.setMappingQuality(mappingQualityCodec.read(bis));

			if (storeMappedQualityScores) {
				boolean hasQS = bis.readBit();
				if (hasQS) {
					bis.alignToByte();
					bisByteAligned = true;
					if (hasQS) {
						byte[] scores = new byte[readLen];
						bis.readAlignedBytes(scores);
						record.setQualityScores(scores);
					}
				}
			}

		} else {
			boolean hasQS = bis.readBit();
			bis.alignToByte();
			bisByteAligned = true;

			byte[] bases = new byte[readLen];
			bis.readAlignedBytes(bases);
			record.setReadBases(bases);

			if (hasQS) {
				byte[] scores = new byte[readLen];
				bis.readAlignedBytes(scores);
				record.setQualityScores(scores);
			}
		}

		timeTick();

		if (!deferredTags.isEmpty()) {
			bis.alignToByte();

			for (int i = 0; i < deferredTags.size(); i++) {
				String tagKeyAndType = deferredTags.get(i);
				int len = tagLens.get(i);
				byte[] value = new byte[len];
				bis.readAlignedBytes(value);

				ReadTag tag = ReadTag.deriveTypeFromKeyAndType(tagKeyAndType,
						ReadTag.restoreValueFromByteArray(tagKeyAndType.charAt(3), value));
				record.tags.add(tag);
			}
		}

		timeTick();

		if (recordCounter % 99999 == 0) {
			StringBuffer sb = new StringBuffer("CramRecordCodec read time: ") ;
			for (int i = 0; i < timeIndex; i++) {
				sb.append(String.format("\t%.3f", (timePoints[i]) / 1000000f));
			}
			log.info(sb.toString()) ;
//			log.info("Total read time: "+readTime);
			Arrays.fill(timePoints, 0);
			readTime = 0 ;
		}

		// marker = bis.readLongBits(debugRecordEndMarkerLen);
		// if (marker != debugRecordEndMarker) {
		// System.out.println(record.toString());
		// throw new
		// RuntimeException("Debug marker for end of record not found.");
		// }
		
		readTime += System.currentTimeMillis()-startMillis ;
		return record;
	}

	private void checkMarker(BitInputStream bis, CramRecord record) throws IOException {
		long marker = bis.readLongBits(debugRecordEndMarkerLen);
		if (marker != debugRecordEndMarker) {
			System.err.println("Record counter=" + recordCounter);
			System.err.println("Record so far: " + record.toString());
			throw new RuntimeException("Oops.");
		}

	}

	private long recordCounter = 0;
	private List<byte[]> deferredByteArrays = new ArrayList<byte[]>();

	@Override
	public long write(BitOutputStream bos, CramRecord record) throws IOException {
		recordCounter++;

		long len = 0L;

		len += flagsCodec.write(bos, record.getFlags());

		if (preserveReadNames)
			len += readNameCodec.write(bos, record.getReadName().getBytes());

		if (!record.isLastFragment()) {
			if (record.getRecordsToNextFragment() > 0) {
				len += recordsToNextFragmentCodec.write(bos, record.getRecordsToNextFragment());
			} else {

				CramRecord mate = record.next == null ? record.previous : record.next;
				bos.write(mate.isReadMapped());
				bos.write(mate.isNegativeStrand());
				bos.write(mate.isFirstInPair());
				if (!preserveReadNames) {
					len += readNameCodec.write(bos, record.getReadName().getBytes());
				}
				len += writeZeroTerminatedString(mate.getSequenceName(), heapByteCodec, bos);
				len += writeZeroTerminatedString(String.valueOf(mate.getAlignmentStart()), heapByteCodec, bos);
				len += writeZeroTerminatedString(String.valueOf(record.insertSize), heapByteCodec, bos);
			}
		}

		if (record.getReadLength() != defaultReadLength) {
			bos.write(true);
			len += readlengthCodec.write(bos, record.getReadLength());
		} else
			bos.write(false);
		len++;

		len += readGroupCodec.write(bos, record.getReadGroupID());

		len += inSeqPosCodec.write(bos, record.getAlignmentStart() - prevPosInSeq);
		prevPosInSeq = record.getAlignmentStart();

		deferredByteArrays.clear();
		if (record.tags == null || record.tags.isEmpty())
			len += tagCountCodec.write(bos, (byte) 0);
		else {
			len += tagCountCodec.write(bos, (byte) record.tags.size());
			for (ReadTag tag : record.tags) {
				len += tagKeyAndTypeCodec.write(bos, tag.getKeyAndType());

				BitCodec<byte[]> tagValueCodec = tagCodecMap.get(tag.getKeyAndType());
				if (tagValueCodec != null) {
					len += tagValueCodec.write(bos, tag.getValueAsByteArray());
				} else {
					BitCodec<Integer> lenCodec = tagValueByteLenCodec.get(tag.getKeyAndType());
					byte[] value = tag.getValueAsByteArray();
					len += lenCodec.write(bos, value.length);

					deferredByteArrays.add(value);
				}
			}
		}

		boolean bosByteAligned = false;

		if (record.isReadMapped()) {
			if (record.getAlignmentStart() - prevPosInSeq < 0) {
				log.error("Negative relative position in sequence: prev=" + prevPosInSeq);
				log.error(record.toString());
			}
			List<ReadFeature> vars = record.getReadFeatures();
			if (vars == null || vars.isEmpty())
				bos.write(false);
			else {
				bos.write(true);
				len += variationsCodec.write(bos, vars);
			}
			len++;

			len += mappingQualityCodec.write(bos, record.getMappingQuality());

			if (storeMappedQualityScores) {
				if (record.getQualityScores() == null || record.getQualityScores().length == 0) {
					bos.write(false);
				} else {
					bos.write(true);
				}
				len++;

				len += bos.alignToByte();
				bosByteAligned = true;

				if (record.getQualityScores() != null && record.getQualityScores().length != 0) {
					bos.write(record.getQualityScores());
					len += 8 * record.getQualityScores().length;
				}
			}

		} else {
			if (record.getAlignmentStart() - prevPosInSeq < 0) {
				log.error("Negative relative position in sequence: prev=" + prevPosInSeq);
				log.error(record.toString());
			}

			if (record.getQualityScores() == null || record.getQualityScores().length == 0) {
				bos.write(false);
			} else {
				bos.write(true);
			}
			len++;

			len += bos.alignToByte();
			bosByteAligned = true;

			bos.write(record.getReadBases());
			len += 8 * record.getReadBases().length;

			if (record.getQualityScores() != null && record.getQualityScores().length != 0) {
				bos.write(record.getQualityScores());
				len += 8 * record.getQualityScores().length;
			}
		}
		if (deferredByteArrays != null && !deferredByteArrays.isEmpty()) {
			len += bos.alignToByte();
			bosByteAligned = true;

			for (byte[] bytes : deferredByteArrays) {
				bos.write(bytes);
				len += 8 * bytes.length;
			}
		}

		// bos.write(debugRecordEndMarker, debugRecordEndMarkerLen);

		return len;
	}

	@Override
	public long numberOfBits(CramRecord record) {
		try {
			return write(NullBitOutputStream.INSTANCE, record);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	private static long writeZeroTerminatedString(String string, BitCodec<Byte> codec, BitOutputStream bos)
			throws IOException {
		long len = 0;
		for (byte b : string.getBytes()) {
			len += codec.write(bos, b);
		}

		len += codec.write(bos, (byte) 0);
		return len;
	}

	private static final int maxBufferSize = 1024;
	private static java.nio.ByteBuffer byteBuffer = java.nio.ByteBuffer.allocate(maxBufferSize);

	private static String readZeroTerminatedString(BitCodec<Byte> codec, BitInputStream bis) throws IOException {
		byteBuffer.clear();
		for (int i = 0; i < maxBufferSize; i++) {
			byte b = codec.read(bis);
			if (b == 0)
				break;
			byteBuffer.put(b);
		}
		if (byteBuffer.position() >= maxBufferSize)
			throw new RuntimeException("Buffer overflow while reading string. ");

		byteBuffer.flip();
		byte[] bytes = new byte[byteBuffer.limit()];
		byteBuffer.get(bytes);
		return new String(bytes);
	}
}
