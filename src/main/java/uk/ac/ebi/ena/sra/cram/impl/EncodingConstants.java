package uk.ac.ebi.ena.sra.cram.impl;

import uk.ac.ebi.ena.sra.cram.Utils;

public class EncodingConstants {
	public static final int INTRA_READ_GOLOMB_M = 8;
	public static final int INTRA_READ_GOLOMB_LOG2M = Utils.mostSignificantBit(INTRA_READ_GOLOMB_M);
	public static final int INTER_READ_GOLOMB_M = 2;
	public static final int INTER_READ_GOLOMB_LOG2M = Utils.mostSignificantBit(INTER_READ_GOLOMB_M);
	public static final int DELETION_LENGTH_GOLOMB_M = 4;
	public static final int DELETION_LENGTH_GOLOMB_LOG2M = Utils.mostSignificantBit(DELETION_LENGTH_GOLOMB_M);
}
