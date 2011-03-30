package uk.ac.ebi.ena.sra.cram.impl;

import java.io.IOException;
import java.util.List;

import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.encoding.BaseChange;
import uk.ac.ebi.ena.sra.cram.encoding.BaseChangeCodec;
import uk.ac.ebi.ena.sra.cram.encoding.BaseSequenceCodec;
import uk.ac.ebi.ena.sra.cram.encoding.BitCodec;
import uk.ac.ebi.ena.sra.cram.encoding.GolombRiceCodec;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;
import uk.ac.ebi.ena.sra.cram.format.Variation;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class CramRecordWriter {

	private BitOutputStream bitOutputStream;
	private CramRecord prevRecord;
	private int prevInReadPos = 0;
	private BitCodec<Long> posInSeqCodec = new GolombRiceCodec(1 << 2);
	private BitCodec<Long> posInReadCodec = new GolombRiceCodec(1 << 3);
	private BitCodec<Long> delLenCodec = new GolombRiceCodec(1 << 3);
	private BitCodec<BaseChange> baseChangeCodec = new BaseChangeCodec();
	private BitCodec<byte[]> basesCodec = new BaseSequenceCodec(
			BaseSequenceCodec.BaseCodecType.FLAT, "SACGTN".getBytes());

	public CramRecordWriter(BitOutputStream os) {
		bitOutputStream = os;
	}

	public void writeCramRecord(CramRecord record) throws IOException {
		if (prevRecord == null)
			writeAlignmentStart(record.getAlignmentStart());
		else
			writeAlignmentStart(record.getAlignmentStart()
					- prevRecord.getAlignmentStart());
		prevRecord = record;
		prevInReadPos = 0;

		appendBit(record.isNegativeStrand());
		if (!record.isPerfectMatch()) {
			appendBit(true);

			List<Variation> vars = Utils.sortVariationsByPosition(record);
			for (Variation v : vars) {
				switch (v.getOperator()) {
				case 'I':
					appendInsertion((InsertionVariation) v);
					break;
				case 'S':
					appendSubstitution((SubstitutionVariation) v);
					break;
				case 'D':
					appendDeletion((DeletionVariation) v);
					break;

				default:
					break;
				}
			}
			appendBit(true);
			appendBit(false);
		} else
			appendBit(false);

	}

	public void flush() throws IOException {
		bitOutputStream.flush();
	}

	private void writeAlignmentStart(long start) throws IOException {
		posInSeqCodec.write(bitOutputStream, (long) start);
	}

	private void appendBit(boolean bit) throws IOException {
		bitOutputStream.writeBits(bit ? 1 : 0, 1);
	}

	private void appendSubstitution(SubstitutionVariation variation)
			throws IOException {
		bitOutputStream.writeBits(0, 1);
		if (variation.getPosition() == 0)
			throw new IllegalArgumentException("Variation position is zero.");

		posInReadCodec.write(bitOutputStream,
				(long) (variation.getPosition() - prevInReadPos));
		baseChangeCodec
				.write(bitOutputStream,
						new BaseChange(variation.getRefernceBase(), variation
								.getBase()));

		prevInReadPos = variation.getPosition();
	}

	private void appendInsertion(InsertionVariation variation)
			throws IOException {
		bitOutputStream.writeBits(1, 1);
		bitOutputStream.writeBits(1, 1);
		bitOutputStream.writeBits(0, 1);
		if (variation.getPosition() == 0)
			throw new IllegalArgumentException("Variation position is zero.");

		posInReadCodec.write(bitOutputStream,
				(long) (variation.getPosition() - prevInReadPos));

		basesCodec.write(bitOutputStream, variation.getSequence());
		prevInReadPos = variation.getPosition();
	}

	private void appendDeletion(DeletionVariation variation) throws IOException {
		bitOutputStream.writeBits(1, 1);
		bitOutputStream.writeBits(1, 1);
		bitOutputStream.writeBits(1, 1);
		if (variation.getPosition() == 0)
			throw new IllegalArgumentException("Variation position is zero.");

		posInReadCodec.write(bitOutputStream,
				(long) (variation.getPosition() - prevInReadPos));
		delLenCodec.write(bitOutputStream, (long) variation.getLength());
		prevInReadPos = variation.getPosition();
	}
}
