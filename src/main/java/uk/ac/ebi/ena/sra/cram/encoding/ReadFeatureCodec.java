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
import java.util.List;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.format.BaseQualityScore;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.InsertBase;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.ReadBase;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.NullBitOutputStream;

public class ReadFeatureCodec implements BitCodec<List<ReadFeature>> {
	public BitCodec<Long> inReadPosCodec;
	public BitCodec<ReadBase> readBaseCodec;
	public BitCodec<SubstitutionVariation> substitutionCodec;
	public BitCodec<InsertionVariation> insertionCodec;
	public BitCodec<DeletionVariation> deletionCodec;
	public BitCodec<BaseQualityScore> baseQSCodec;
	public BitCodec<InsertBase> insertBaseCodec;

	public BitCodec<Byte> featureOperationCodec;

	private static Logger log = Logger.getLogger(ReadFeatureCodec.class);

	@Override
	public List<ReadFeature> read(BitInputStream bis) throws IOException {
		List<ReadFeature> list = new ArrayList<ReadFeature>();
		byte op;
		int prevPos = 0;
		while ((op = featureOperationCodec.read(bis)) != ReadFeature.STOP_OPERATOR) {
			ReadFeature feature = null;
			int pos = prevPos + inReadPosCodec.read(bis).intValue();
			prevPos = pos;
			switch (op) {
			case ReadBase.operator:
				ReadBase readBase = readBaseCodec.read(bis);
				readBase.setPosition(pos);
				feature = readBase;
				break;
			case SubstitutionVariation.operator:
				SubstitutionVariation sub = substitutionCodec.read(bis);
				sub.setPosition(pos);
				feature = sub;
				break;
			case InsertionVariation.operator:
				InsertionVariation ins = insertionCodec.read(bis);
				ins.setPosition(pos);
				feature = ins;
				break;
			case DeletionVariation.operator:
				DeletionVariation del = deletionCodec.read(bis);
				del.setPosition(pos);
				feature = del;
				break;
			case InsertBase.operator:
				InsertBase ib = insertBaseCodec.read(bis);
				ib.setPosition(pos);
				feature = ib;
				break;
			case BaseQualityScore.operator:
				BaseQualityScore bqs = baseQSCodec.read(bis);
				bqs.setPosition(pos);
				feature = bqs;
				break;

			default:
				throw new RuntimeException("Unknown read feature operator: "
						+ (char) op);
			}
			list.add(feature);
		}

		return list;
	}

	@Override
	public long write(BitOutputStream bos, List<ReadFeature> features)
			throws IOException {

		long len = 0L;
		int prevPos = 0;
		for (ReadFeature feature : features) {
			len += featureOperationCodec.write(bos, feature.getOperator());

			len += inReadPosCodec.write(bos, (long) feature.getPosition()
					- prevPos);

			prevPos = feature.getPosition();
			switch (feature.getOperator()) {
			case ReadBase.operator:
				len += readBaseCodec.write(bos, (ReadBase) feature);
				break;
			case SubstitutionVariation.operator:

				len += substitutionCodec.write(bos,
						(SubstitutionVariation) feature);
				break;
			case InsertionVariation.operator:
				len += insertionCodec.write(bos, (InsertionVariation) feature);
				break;
			case DeletionVariation.operator:
				len += deletionCodec.write(bos, (DeletionVariation) feature);
				break;
			case InsertBase.operator:
				len += insertBaseCodec.write(bos, (InsertBase) feature);
				break;
			case BaseQualityScore.operator:
				len += baseQSCodec.write(bos, (BaseQualityScore) feature);
				break;

			default:
				throw new RuntimeException("Unknown read feature operator: "
						+ (char) feature.getOperator());
			}
		}
		len += featureOperationCodec.write(bos, ReadFeature.STOP_OPERATOR);

		return len;
	}

	@Override
	public long numberOfBits(List<ReadFeature> features) {
		try {
			return write(NullBitOutputStream.INSTANCE, features);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
