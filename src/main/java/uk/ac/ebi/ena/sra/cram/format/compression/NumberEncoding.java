package uk.ac.ebi.ena.sra.cram.format.compression;

import java.util.TreeMap;

public enum NumberEncoding {
	UNARY, BETA, GAMMA, GOLOMB, GOLOMB_RICE, SUBEXP;
	
	private static TreeMap<String, NumberEncoding> tagmap = new TreeMap<String, NumberEncoding>() ;
	static {
		tagmap.put("UN", UNARY) ;
	}
}
