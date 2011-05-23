package uk.ac.ebi.ena.sra.compression.huffman;

public class HuffmanLeaf<T> extends HuffmanTree<T> {
	// user object, attached to the leaf: 
	public final T value; 

	public HuffmanLeaf(int freq, T val) {
		super(freq);
		value = val;
	}
}
