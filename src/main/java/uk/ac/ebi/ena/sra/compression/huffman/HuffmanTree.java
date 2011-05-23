package uk.ac.ebi.ena.sra.compression.huffman;

public abstract class HuffmanTree<T> implements Comparable<HuffmanTree<T>> {
	public final int frequency; 

	public HuffmanTree(int freq) {
		frequency = freq;
	}

	public int compareTo(HuffmanTree<T> tree) {
		return frequency - tree.frequency;
	}
}
