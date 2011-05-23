package uk.ac.ebi.ena.sra.compression.huffman;

import java.io.PrintStream;
import java.util.PriorityQueue;

public class HuffmanCode {

	public static <T> HuffmanTree<T> buildTree(int[] charFreqs, T[] values) {
		PriorityQueue<HuffmanTree<T>> queue = new PriorityQueue<HuffmanTree<T>>();

		for (int i = 0; i < charFreqs.length; i++)
			if (charFreqs[i] > 0)
				queue.offer(new HuffmanLeaf<T>(charFreqs[i], values[i]));

		while (queue.size() > 1) {
			HuffmanTree<T> a = queue.poll();
			HuffmanTree<T> b = queue.poll();

			queue.offer(new HuffmanNode<T>(a, b));
		}
		return queue.poll();
	}

	public static void printTree(HuffmanTree<?> tree, StringBuffer prefix,
			PrintStream ps) {
		if (tree instanceof HuffmanLeaf) {
			HuffmanLeaf<?> leaf = (HuffmanLeaf<?>) tree;

			ps.println(leaf.value + "\t" + leaf.frequency + "\t" + prefix);

		} else if (tree instanceof HuffmanNode) {
			HuffmanNode<?> node = (HuffmanNode<?>) tree;

			// traverse left
			prefix.append('0');
			printTree(node.left, prefix, ps);
			prefix.deleteCharAt(prefix.length() - 1);

			// traverse right
			prefix.append('1');
			printTree(node.right, prefix, ps);
			prefix.deleteCharAt(prefix.length() - 1);
		}
	}
}
