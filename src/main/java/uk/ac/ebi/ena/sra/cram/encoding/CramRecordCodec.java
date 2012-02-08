package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;
import java.util.List;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.ReadAnnotation;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.NullBitOutputStream;

public class CramRecordCodec implements BitCodec<CramRecord> {
	public BitCodec<Long> inSeqPosCodec;
	public BitCodec<Long> recordsToNextFragmentCodec;
	public BitCodec<Long> readlengthCodec;
	public BitCodec<List<ReadFeature>> variationsCodec;
	public SequenceBaseProvider sequenceBaseProvider;
	// public BitCodec<byte[]> basesCodec;
	// public BitCodec<byte[]> qualitiesCodec;

	public BitCodec<Byte> baseCodec;
	public BitCodec<Byte> qualityCodec;

	public String sequenceName;
	public long prevPosInSeq = 1L;
	public long defaultReadLength = 0L;

	public BitCodec<ReadAnnotation> readAnnoCodec;
	public BitCodec<Integer> readGroupCodec;
	public BitCodec<Long> nextFragmentIDCodec;

	public BitCodec<Byte> mappingQualityCodec;

	public boolean storeMappedQualityScores = false;
	public BitCodec<Byte> heapByteCodec;

	private static Logger log = Logger.getLogger(CramRecordCodec.class);

	@Override
	public CramRecord read(BitInputStream bis) throws IOException {
		CramRecord record = new CramRecord();

		record.setFirstInPair(bis.readBit());
		record.setProperPair(bis.readBit());
		record.setDuplicate(bis.readBit());

		record.setNegativeStrand(bis.readBit());
		record.setLastFragment(bis.readBit());
		if (!record.isLastFragment()) {
			if (bis.readBit()) {
				record.setRecordsToNextFragment(recordsToNextFragmentCodec.read(bis));
			} else {
				CramRecord mate = new CramRecord();
				mate.setReadMapped(bis.readBit());
				mate.setNegativeStrand(bis.readBit());
				mate.setFirstInPair(bis.readBit());
				mate.setReadName(readZeroTerminatedString(heapByteCodec, bis));
				mate.setSequenceName(readZeroTerminatedString(heapByteCodec, bis));
				mate.setAlignmentStart(Long.valueOf(readZeroTerminatedString(heapByteCodec, bis)));
				mate.setFirstInPair(!record.isFirstInPair());
				if (record.isFirstInPair())
					record.next = mate;
				else
					record.previous = mate;

				record.setReadName(mate.getReadName());
			}
		}

		int readLen;
		if (bis.readBit())
			readLen = readlengthCodec.read(bis).intValue();
		else
			readLen = (int) defaultReadLength;
		record.setReadLength(readLen);

		boolean readMapped = bis.readBit();
		record.setReadMapped(readMapped);
		if (readMapped) {
			long position = prevPosInSeq + inSeqPosCodec.read(bis);
			prevPosInSeq = position;
			record.setAlignmentStart(position);

			boolean imperfectMatch = bis.readBit();
			record.setPerfectMatch(!imperfectMatch);
			if (imperfectMatch) {
				List<ReadFeature> features = variationsCodec.read(bis);
				record.setReadFeatures(features);
			}

			if (storeMappedQualityScores) {
				byte[] scores = new byte[readLen];
				readNonEmptyByteArray(bis, scores, qualityCodec);
				record.setQualityScores(scores);
			}

			record.setMappingQuality(mappingQualityCodec.read(bis));
		} else {
			long position = prevPosInSeq + inSeqPosCodec.read(bis);
			prevPosInSeq = position;
			record.setAlignmentStart(position);

			byte[] bases = new byte[readLen];
			readNonEmptyByteArray(bis, bases, baseCodec);
			record.setReadBases(bases);

			byte[] scores = new byte[readLen];
			readNonEmptyByteArray(bis, scores, qualityCodec);
			record.setQualityScores(scores);
		}

		record.setReadGroupID(readGroupCodec.read(bis));

		// if (bis.readBit()) {
		// List<ReadAnnotation> anns = new ArrayList<ReadAnnotation>();
		// do {
		// anns.add(readAnnoCodec.read(bis));
		// } while (bis.readBit());
		// if (!anns.isEmpty())
		// record.setAnnotations(anns);
		// }

		return record;
	}

	@Override
	public long write(BitOutputStream bos, CramRecord record) throws IOException {
		long len = 0L;

		bos.write(record.isFirstInPair());
		len++;
		bos.write(record.isProperPair());
		len++;
		bos.write(record.isDuplicate());
		len++;

		bos.write(record.isNegativeStrand());
		len++;

		bos.write(record.isLastFragment());
		len++;

		if (!record.isLastFragment()) {
			if (record.getRecordsToNextFragment() > 0) {
				bos.write(true);
				len += recordsToNextFragmentCodec.write(bos, record.getRecordsToNextFragment());
			} else {
				bos.write(false);

				CramRecord mate = record.next == null ? record.previous : record.next;
				bos.write(mate.isReadMapped());
				bos.write(mate.isNegativeStrand());
				bos.write(mate.isFirstInPair());
				len += writeZeroTerminatedString(record.getReadName(), heapByteCodec, bos);
				len += writeZeroTerminatedString(mate.getSequenceName(), heapByteCodec, bos);
				len += writeZeroTerminatedString(String.valueOf(mate.getAlignmentStart()), heapByteCodec, bos);
			}
			len++;
		}

		if (record.getReadLength() != defaultReadLength) {
			bos.write(true);
			len += readlengthCodec.write(bos, record.getReadLength());
		} else
			bos.write(false);
		len++;

		if (record.isReadMapped()) {
			bos.write(true);
			len++;

			if (record.getAlignmentStart() - prevPosInSeq < 0) {
				log.error("Negative relative position in sequence: prev=" + prevPosInSeq);
				log.error(record.toString());
			}
			len += inSeqPosCodec.write(bos, record.getAlignmentStart() - prevPosInSeq);
			prevPosInSeq = record.getAlignmentStart();

			List<ReadFeature> vars = record.getReadFeatures();
			if (vars == null || vars.isEmpty())
				bos.write(false);
			else {
				bos.write(true);
				len += variationsCodec.write(bos, vars);
			}
			len++;

			if (storeMappedQualityScores)
				len += writeNonEmptyByteArray(bos, record.getQualityScores(), qualityCodec);

			mappingQualityCodec.write(bos, record.getMappingQuality());
		} else {
			bos.write(false);
			len++;

			if (record.getAlignmentStart() - prevPosInSeq < 0) {
				log.error("Negative relative position in sequence: prev=" + prevPosInSeq);
				log.error(record.toString());
			}
			len += inSeqPosCodec.write(bos, record.getAlignmentStart() - prevPosInSeq);
			prevPosInSeq = record.getAlignmentStart();

			len += writeNonEmptyByteArray(bos, record.getReadBases(), baseCodec);
			len += writeNonEmptyByteArray(bos, record.getQualityScores(), qualityCodec);
		}

		len += readGroupCodec.write(bos, record.getReadGroupID());

		// Collection<ReadAnnotation> annotations = record.getAnnotations();
		// if (annotations == null || annotations.isEmpty()) {
		// bos.write(false);
		// len++;
		// } else {
		// for (ReadAnnotation a : annotations) {
		// bos.write(true);
		// len++;
		// len += readAnnoCodec.write(bos, a);
		// }
		// bos.write(false);
		// len++;
		// }

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

	private static int writeNonEmptyByteArray(BitOutputStream bos, byte[] array, BitCodec<Byte> codec)
			throws IOException {
		if (array == null || array.length == 0)
			throw new RuntimeException("Expecting a non-empty array.");

		int len = 0;
		for (byte b : array)
			len += codec.write(bos, b);
		return len;
	}

	private static byte[] readNonEmptyByteArray(BitInputStream bis, byte[] array, BitCodec<Byte> codec)
			throws IOException {
		for (int i = 0; i < array.length; i++)
			array[i] = codec.read(bis);

		return array;
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
