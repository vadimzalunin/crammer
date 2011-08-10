package uk.ac.ebi.ena.sra.cram.encoding;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.cram.format.ReadAnnotation;

public class ReadAnnotationCodec extends HuffmanCodec<ReadAnnotation> {

	public ReadAnnotationCodec(ReadAnnotation[] anns, int[] freqs) {
		super(HuffmanCode.buildTree(freqs, anns));
	}

}
