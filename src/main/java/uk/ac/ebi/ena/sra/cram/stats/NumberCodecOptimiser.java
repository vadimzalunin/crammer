package uk.ac.ebi.ena.sra.cram.stats;

import java.util.ArrayList;
import java.util.List;

import uk.ac.ebi.ena.sra.cram.format.compression.CramCompressionException;
import uk.ac.ebi.ena.sra.cram.format.compression.EncodingAlgorithm;
import uk.ac.ebi.ena.sra.cram.format.compression.NumberCodecFactory;
import uk.ac.ebi.ena.sra.cram.format.compression.NumberCodecStub;

public class NumberCodecOptimiser {

	private NumberCodecStub[] stubs;
	private long lengths[];

	public NumberCodecOptimiser() throws CramCompressionException {
		List<NumberCodecStub> list = new ArrayList<NumberCodecStub>();
		NumberCodecStub stub = null;
		for (int i = 2; i < 20; i++) {
			stub = NumberCodecFactory.createStub(EncodingAlgorithm.GOLOMB);
			stub.initFromString(i + ",0,1");
			list.add(stub);
		}

		for (int i = 1; i < 20; i++) {
			stub = NumberCodecFactory.createStub(EncodingAlgorithm.GOLOMB_RICE);
			stub.initFromString(i + ",0,1");
			list.add(stub);
		}

		stub = NumberCodecFactory.createStub(EncodingAlgorithm.GAMMA);
		stub.initFromString("1,0");
		list.add(stub);

		for (int i = 1; i < 20; i++) {
			stub = NumberCodecFactory.createStub(EncodingAlgorithm.SUBEXP);
			stub.initFromString(i + ",0,1");
			list.add(stub);
		}

		stub = NumberCodecFactory.createStub(EncodingAlgorithm.UNARY);
		stub.initFromString("0,1");
		list.add(stub);

		stubs = (NumberCodecStub[]) list.toArray(new NumberCodecStub[list
				.size()]);
		lengths = new long[stubs.length];
	}

	public void addValue(long value, long count) {
		for (int i = 0; i < stubs.length; i++) {
			NumberCodecStub stub = stubs[i];
			lengths[i] += stub.numberOfBits(value) * count;
		}
	}

	public NumberCodecStub getMinLengthStub() {
		int minLengtIndex = 0;
		for (int i = 1; i < stubs.length; i++) {
			if (lengths[minLengtIndex] > lengths[i])
				minLengtIndex = i;
		}

		return stubs[minLengtIndex];
	}

	public long getMinLength() {
		int minLengtIndex = 0;
		for (int i = 1; i < stubs.length; i++) {
			if (lengths[minLengtIndex] > lengths[i])
				minLengtIndex = i;
		}

		return lengths[minLengtIndex];
	}

	@Override
	public String toString() {
		// TODO Auto-generated method stub
		return super.toString();
	}

}
