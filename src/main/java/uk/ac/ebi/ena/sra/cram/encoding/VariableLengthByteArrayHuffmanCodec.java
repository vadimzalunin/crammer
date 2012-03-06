package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class VariableLengthByteArrayHuffmanCodec implements BitCodec<byte[]> {
	private HuffmanCodec<Byte> byteCodec;
	private HuffmanCodec<Integer> lenCodec;

	public VariableLengthByteArrayHuffmanCodec(byte[] alphabet, int[] freqs, int[] lenAlphabet, int[] lenFreqs) {

		int maxLen = 0;
		for (int len : lenAlphabet)
			if (maxLen < len)
				maxLen = len;

		HuffmanTree<Byte> tree = HuffmanCode.buildTree(freqs, Utils.autobox(alphabet));
		byteCodec = new HuffmanCodec<Byte>(tree);

		HuffmanTree<Integer> lenTree = HuffmanCode.buildTree(lenFreqs, Utils.autobox(lenAlphabet));
		lenCodec = new HuffmanCodec<Integer>(lenTree);
	}

	@Override
	public byte[] read(BitInputStream bis) throws IOException {
		Integer len = lenCodec.read(bis) ;
		byte[] sequence = new byte[len];

		for (int i = 0; i < sequence.length; i++)
			sequence[i] = byteCodec.read(bis);

		return sequence;
	}

	@Override
	public long write(BitOutputStream bos, byte[] bytes) throws IOException {
		long len = 0;

		len += lenCodec.write(bos, bytes.length);

		for (byte b : bytes)
			len += byteCodec.write(bos, b);

		return len;
	}

	@Override
	public long numberOfBits(byte[] bytes) {
		long len = 0;

		len += lenCodec.numberOfBits(bytes.length);

		for (byte b : bytes)
			len += byteCodec.numberOfBits(b);

		return len;
	}

}
