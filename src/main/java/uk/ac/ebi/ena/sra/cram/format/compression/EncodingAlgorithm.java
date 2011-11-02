package uk.ac.ebi.ena.sra.cram.format.compression;

import java.io.Serializable;

public enum EncodingAlgorithm implements Serializable{
	NULL, UNARY, BETA, GAMMA, GOLOMB, GOLOMB_RICE, SUBEXP;
}
