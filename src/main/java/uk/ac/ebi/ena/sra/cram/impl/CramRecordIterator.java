package uk.ac.ebi.ena.sra.cram.impl;

import java.io.EOFException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.encoding.BaseChange;
import uk.ac.ebi.ena.sra.cram.encoding.BaseChangeCodec;
import uk.ac.ebi.ena.sra.cram.encoding.BaseSequenceCodec;
import uk.ac.ebi.ena.sra.cram.encoding.BitCodec;
import uk.ac.ebi.ena.sra.cram.encoding.GolombRiceCodec;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;

public class CramRecordIterator implements Iterator<CramRecord> {
	private Long nextPos = 1L;
	private BitInputStream decodingStream;
	private SequenceBaseProvider sequenceBaseProvider;
	private CramRecord prevRecord;
	private int prevInReadPos;
	private String seqName;

	private BitCodec<Long> posInSeqCodec = new GolombRiceCodec(1 << 2);
	private BitCodec<Long> posInReadCodec = new GolombRiceCodec(1 << 3);
	private BitCodec<Long> delLenCodec = new GolombRiceCodec(1 << 3);
	private BitCodec<BaseChange> baseChangeCodec = new BaseChangeCodec();
	private BitCodec<byte[]> basesCodec = new BaseSequenceCodec(
			BaseSequenceCodec.BaseCodecType.FLAT, "SACGTN".getBytes());

	public CramRecordIterator(BitInputStream is,
			SequenceBaseProvider sequenceBaseProvider) {
		super();
		decodingStream = is;
		this.sequenceBaseProvider = sequenceBaseProvider;
	}

	public void setSequenceName(String seqName) {
		this.seqName = seqName;
	}

	@Override
	public boolean hasNext() {
		try {
			nextPos = posInSeqCodec.read(decodingStream);
		} catch (EOFException e) {
			return false;
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		return true;
	}

	@Override
	public CramRecord next() {
		if (nextPos < 0)
			throw new RuntimeException("Next pos is non-positive: " + nextPos);
		CramRecord record = new CramRecord();
		if (prevRecord == null)
			record.setAlignmentStart(nextPos);
		else
			record.setAlignmentStart(prevRecord.getAlignmentStart() + nextPos);
		prevRecord = record;

		nextPos = 0L;

		prevInReadPos = 0;

		try {
			record.setNegativeStrand(decodingStream.readBit());
			record.setPerfectMatch(!decodingStream.readBit());

			if (!record.isPerfectMatch()) {
				long pos = 0;

				do {
					if (!decodingStream.readBit()) {
						// subs:

						Collection<SubstitutionVariation> subs = record
								.getSubstitutionVariations();
						if (record.getSubstitutionVariations() == null) {
							subs = new ArrayList<SubstitutionVariation>();
							record.setSubstitutionVariations(subs);
						}

						while ((pos = posInReadCodec.read(decodingStream)) != 0) {

							SubstitutionVariation v = new SubstitutionVariation();
							prevInReadPos += pos;
							v.setPosition(prevInReadPos);

							BaseChange baseChange = baseChangeCodec
									.read(decodingStream);

							byte refBase = sequenceBaseProvider.getBaseAt(
									seqName,
									record.getAlignmentStart()
											+ v.getPosition() - 2);

							byte base = baseChange.getBaseForReference(refBase);

							v.setRefernceBase(refBase);
							v.setBase(base);
							subs.add(v);
							break;
						}

					} else if (!decodingStream.readBit()) {// stop
						break;
					} else {

						if (!decodingStream.readBit()) {
							// ins:
							Collection<InsertionVariation> insertions = record
									.getInsertionVariations();
							if (insertions == null) {
								insertions = new ArrayList<InsertionVariation>();
								record.setInsertionVariations(insertions);
							}
							while ((pos = posInReadCodec.read(decodingStream)) != 0) {
								byte[] bases = basesCodec.read(decodingStream);
								InsertionVariation v = new InsertionVariation();
								prevInReadPos += pos;
								v.setPosition(prevInReadPos);
								v.setSequence(bases);
								insertions.add(v);
								break;
							}
						} else {
							// dels:
							Collection<DeletionVariation> deletions = record
									.getDeletionVariations();
							if (deletions == null) {
								deletions = new ArrayList<DeletionVariation>();
								record.setDeletionVariations(deletions);
							}
							while ((pos = posInReadCodec.read(decodingStream)) != 0) {
								long length = delLenCodec.read(decodingStream);
								DeletionVariation v = new DeletionVariation();
								prevInReadPos += pos;
								v.setPosition(prevInReadPos);
								v.setLength((int) length);
								deletions.add(v);
								break;
							}
						}
					}
				} while (true);

				// ugly hack:
				if (record.getSubstitutionVariations() != null
						&& record.getSubstitutionVariations().isEmpty())
					record.setSubstitutionVariations(null);

				if (record.getInsertionVariations() != null
						&& record.getInsertionVariations().isEmpty())
					record.setInsertionVariations(null);

				if (record.getDeletionVariations() != null
						&& record.getDeletionVariations().isEmpty())
					record.setDeletionVariations(null);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return record;
	}

	@Override
	public void remove() {
		throw new RuntimeException(
				"Remove method is not supported in cram reader.");
	}

}
