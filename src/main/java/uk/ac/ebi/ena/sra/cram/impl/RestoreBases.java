package uk.ac.ebi.ena.sra.cram.impl;

import java.io.IOException;
import java.util.List;

import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;

public class RestoreBases {

	private SequenceBaseProvider provider;
	private String sequenceName;

	public RestoreBases() {
	}

	public RestoreBases(SequenceBaseProvider provider, String sequenceName) {
		this.provider = provider;
		this.sequenceName = sequenceName;
	}

	public byte[] restoreReadBases(CramRecord record) throws IOException {
		byte[] bases = new byte[(int) record.getReadLength()];

		int posInRead = 1;
		long posInSeq = record.getAlignmentStart() - 1;
		if (record.getReadFeatures() == null
				|| record.getReadFeatures().isEmpty()) {
			for (posInRead = 1; posInRead <= record.getReadLength(); posInRead++)
				bases[posInRead - 1] = provider.getBaseAt(sequenceName,
						posInSeq++);
			return bases;
		}
		List<ReadFeature> variations = record.getReadFeatures();
		for (ReadFeature v : variations) {
			for (; posInRead < v.getPosition(); posInRead++)
				bases[posInRead - 1] = provider.getBaseAt(sequenceName,
						posInSeq++);

			switch (v.getOperator()) {
			case SubstitutionVariation.operator:
				SubstitutionVariation sv = (SubstitutionVariation) v;
				byte refBase = provider.getBaseAt(sequenceName, posInSeq);
				byte base = sv.getBaseChange().getBaseForReference(refBase);
				sv.setBase(base);
				sv.setRefernceBase(refBase);
				bases[posInRead++ - 1] = sv.getBase();
				posInSeq++;
				break;
			case InsertionVariation.operator:
				InsertionVariation iv = (InsertionVariation) v;
				for (int i = 0; i < iv.getSequence().length; i++)
					bases[posInRead++ - 1] = iv.getSequence()[i];

				break;
			case DeletionVariation.operator:
				DeletionVariation dv = (DeletionVariation) v;
				posInSeq += dv.getLength();
				break;

			default:
				throw new RuntimeException("Uknown variation operator: "
						+ v.getOperator());
			}
		}
		for (; posInRead <= record.getReadLength(); posInRead++)
			bases[posInRead - 1] = provider.getBaseAt(sequenceName, posInSeq++);
		return bases;
	}

	public SequenceBaseProvider getProvider() {
		return provider;
	}

	public void setProvider(SequenceBaseProvider provider) {
		this.provider = provider;
	}

	public String getSequenceName() {
		return sequenceName;
	}

	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}
}
