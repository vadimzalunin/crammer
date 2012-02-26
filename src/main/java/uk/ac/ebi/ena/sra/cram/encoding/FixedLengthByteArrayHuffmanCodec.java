package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;
import java.nio.ByteBuffer;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class FixedLengthByteArrayHuffmanCodec implements BitCodec<byte[]> {
	private HuffmanCodec<Byte> byteCodec;
	private ByteBuffer buf;
	private int len;

	public FixedLengthByteArrayHuffmanCodec(byte[] alphabet, int[] freqs, int len) {
		this.len = len;
		buf = ByteBuffer.allocate(len);
		HuffmanTree<Byte> tree = HuffmanCode.buildTree(freqs, Utils.autobox(alphabet));
		byteCodec = new HuffmanCodec<Byte>(tree);
	}

	@Override
	public byte[] read(BitInputStream bis) throws IOException {
		buf.clear();
		for (int i = 0; i < len; i++)
			buf.put(byteCodec.read(bis));

		byte[] sequence = new byte[buf.position()];
		buf.flip();
		buf.get(sequence);
		return sequence;
	}

	@Override
	public long write(BitOutputStream bos, byte[] bytes) throws IOException {
		if (bytes.length != len)
			throw new RuntimeException("Number of bytes in the value is different from " + len);

		long len = 0;
		for (byte b : bytes)
			len += byteCodec.write(bos, b);

		return len;
	}

	@Override
	public long numberOfBits(byte[] bytes) {
		if (bytes.length != len)
			throw new RuntimeException("Number of bytes in the value is different from " + len);

		long len = 0;
		for (byte b : bytes)
			len += byteCodec.numberOfBits(b);

		return len;
	}

}
