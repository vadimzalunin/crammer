package uk.ac.ebi.ena.sra.cram.format;

import java.io.Serializable;

public class CramCompression implements Serializable {

	int inSeqPosGolombLogM;
	int inSeqPosGolombM;
	int inReadPosGolombLogM;
	int inReadPosGolombM;
	int delLengthGolombLogM;
	int delLengthGolombM;

	public int getInSeqPosGolombLogM() {
		return inSeqPosGolombLogM;
	}

	public int getInSeqPosGolombM() {
		return inSeqPosGolombM;
	}

	public int getInReadPosGolombLogM() {
		return inReadPosGolombLogM;
	}

	public int getInReadPosGolombM() {
		return inReadPosGolombM;
	}

	public int getDelLengthGolombLogM() {
		return delLengthGolombLogM;
	}

	public int getDelLengthGolombM() {
		return delLengthGolombM;
	}

	public void setInSeqPosGolombLogM(int inSeqPosGolombLogM) {
		this.inSeqPosGolombLogM = inSeqPosGolombLogM;
		this.inSeqPosGolombM = 1 << inSeqPosGolombLogM;
	}

	public void setInReadPosGolombLogM(int inReadPosGolombLogM) {
		this.inReadPosGolombLogM = inReadPosGolombLogM;
		this.inReadPosGolombM = 1 << inReadPosGolombLogM;
	}

	public void setDelLengthGolombLogM(int delLengthGolombLogM) {
		this.delLengthGolombLogM = delLengthGolombLogM;
		this.delLengthGolombM = 1 << delLengthGolombLogM;
	}

}
