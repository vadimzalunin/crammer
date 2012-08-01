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

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.InsertBase;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.ReadBase;
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

	private static final long calcRefLength(CramRecord record) {
		if (record.getReadFeatures() == null
				|| record.getReadFeatures().isEmpty())
			return record.getReadLength();
		long len = record.getReadLength();
		for (ReadFeature rf : record.getReadFeatures()) {
			switch (rf.getOperator()) {
			case DeletionVariation.operator:
				len += ((DeletionVariation) rf).getLength();
				break;
			case InsertionVariation.operator:
				len -= ((InsertionVariation) rf).getSequence().length;
				break;
			default:
				break;
			}
		}

		return len;
	}

	public byte[] restoreReadBases(CramRecord record) throws IOException {
		int readLength = (int) record.getReadLength();
		byte[] bases = new byte[readLength];
		// byte[] refBases = new byte[readLength * 2];
		byte[] refBases = new byte[(int) calcRefLength(record)];

		int posInRead = 1;
		long alignmentStart = record.getAlignmentStart() - 1;

		int posInSeq = 0;
		provider.copyBases(sequenceName, alignmentStart, refBases.length,
				refBases);
		if (record.getReadFeatures() == null
				|| record.getReadFeatures().isEmpty()) {
			for (posInRead = 1; posInRead <= readLength; posInRead++)
				bases[posInRead - 1] = refBases[posInSeq++];
			
			record.setReadBases(bases) ;
			return bases;
		}
		List<ReadFeature> variations = record.getReadFeatures();
		for (ReadFeature v : variations) {
			for (; posInRead < v.getPosition(); posInRead++)
				bases[posInRead - 1] = provider.getBaseAt(sequenceName,
						alignmentStart + posInSeq++);

			switch (v.getOperator()) {
			case SubstitutionVariation.operator:
				SubstitutionVariation sv = (SubstitutionVariation) v;
				byte refBase = provider.getBaseAt(sequenceName, alignmentStart
						+ posInSeq);
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
			case InsertBase.operator:
				InsertBase ib = (InsertBase) v;
				bases[posInRead++ - 1] = ib.getBase();
				break;

			// default:
			// throw new RuntimeException("Uknown variation operator: "
			// + v.getOperator());
			}
		}
		for (; posInRead <= readLength; posInRead++)
			bases[posInRead - 1] = provider.getBaseAt(sequenceName,
					alignmentStart + posInSeq++);

		// ReadBase overwrites bases:
		for (ReadFeature v : variations) {
			switch (v.getOperator()) {
			case ReadBase.operator:
				ReadBase rb = (ReadBase) v;
				bases[v.getPosition() - 1] = rb.getBase();
				break;
			default:
				break;
			}
		}
		
		record.setReadBases(bases) ;
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
