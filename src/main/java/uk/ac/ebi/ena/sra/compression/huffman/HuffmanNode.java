package uk.ac.ebi.ena.sra.compression.huffman;

public class HuffmanNode<T> extends HuffmanTree<T> {
	public final HuffmanTree<T> left, right;

	public HuffmanNode(HuffmanTree<T> l, HuffmanTree<T> r) {
		super(l.frequency + r.frequency);
		left = l;
		right = r;
	}
}
