package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;
import java.util.List;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
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
	public BitCodec<byte[]> basesCodec;
	public BitCodec<byte[]> qualitiesCodec;
	public String sequenceName;
	public long prevPosInSeq = 1L;
	public long defaultReadLength = 0L;

	private static Logger log = Logger.getLogger(CramRecordCodec.class);

	@Override
	public CramRecord read(BitInputStream bis) throws IOException {
		CramRecord record = new CramRecord();

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
		} else {
			record.setReadBases(basesCodec.read(bis));
			record.setQualityScores(qualitiesCodec.read(bis));
		}

		record.setNegativeStrand(bis.readBit());
		record.setLastFragment(bis.readBit());
		if (!record.isLastFragment())
			record.setRecordsToNextFragment(recordsToNextFragmentCodec
					.read(bis));

		if (bis.readBit())
			record.setReadLength(readlengthCodec.read(bis));
		else
			record.setReadLength(defaultReadLength);

		return record;
	}

	@Override
	public long write(BitOutputStream bos, CramRecord record)
			throws IOException {
		long len = 0L;
		if (record.isReadMapped()) {
			bos.write(true);
			len++;
			if (record.getAlignmentStart() - prevPosInSeq < 0) {
				log.error("Negative relative position in sequence: prev="
						+ prevPosInSeq);
				log.error(record.toString());
			}
			len += inSeqPosCodec.write(bos, record.getAlignmentStart()
					- prevPosInSeq);

			prevPosInSeq = record.getAlignmentStart();
			if (!record.isPerfectMatch()) {
				bos.write(true);
				List<ReadFeature> vars = record.getReadFeatures();
				len += variationsCodec.write(bos, vars);
			} else
				bos.write(false);
			len++;
		} else {
			bos.write(false);
			len++;

			len += basesCodec.write(bos, record.getReadBases());
			len += qualitiesCodec.write(bos, record.getQualityScores());

			// throw new RuntimeException("Unmapped reads are not supported.");
		}

		bos.write(record.isNegativeStrand());
		len++;
		bos.write(record.isLastFragment());
		len++;
		if (!record.isLastFragment()) {
			len += recordsToNextFragmentCodec.write(bos,
					record.getRecordsToNextFragment());
		}

		if (record.getReadLength() != defaultReadLength) {
			bos.write(true);
			len += readlengthCodec.write(bos, record.getReadLength());
		} else
			bos.write(false);
		len++;

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
}
