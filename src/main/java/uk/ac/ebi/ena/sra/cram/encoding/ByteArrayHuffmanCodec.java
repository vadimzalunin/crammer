package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;
import java.nio.ByteBuffer;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class ByteArrayHuffmanCodec implements BitCodec<byte[]> {
	private HuffmanCodec<Byte> byteCodec;
	private ByteBuffer buf = ByteBuffer.allocate(1000);
	private byte stopByte;

	public ByteArrayHuffmanCodec(byte[] alphabet, int[] freqs, byte stopByte) {
		this.stopByte = stopByte;
		HuffmanTree<Byte> tree = HuffmanCode.buildTree(freqs,
				Utils.autobox(alphabet));
		byteCodec = new HuffmanCodec<Byte>(tree);
	}

	@Override
	public byte[] read(BitInputStream bis) throws IOException {
		buf.clear();
		byte b;
		while ((b = byteCodec.read(bis)) != stopByte)
			buf.put(b);

		byte[] sequence = new byte[buf.position()];
		buf.flip() ;
		buf.get(sequence);
		return sequence;
	}

	@Override
	public long write(BitOutputStream bos, byte[] bytes) throws IOException {
		if (bytes.length == 0)
			return 0;

		long len = 0;
		for (byte b : bytes)
			len += byteCodec.write(bos, b);

		if (bytes[bytes.length - 1] != stopByte)
			len += byteCodec.write(bos, stopByte);

		return len;
	}

	@Override
	public long numberOfBits(byte[] bytes) {
		if (bytes.length == 0)
			return 0;

		long len = 0;
		for (byte b : bytes)
			len += byteCodec.numberOfBits(b);

		if (bytes[bytes.length - 1] != stopByte)
			len += byteCodec.numberOfBits(stopByte);

		return len;
	}

}
