package uk.ac.ebi.ena.sra.cram.format.compression;

import java.io.Serializable;

public enum EncodingAlgorithm implements Serializable{
	UNARY, BETA, GAMMA, GOLOMB, GOLOMB_RICE, SUBEXP, HUFFMAN;
}
