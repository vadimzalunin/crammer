package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;
import java.nio.ByteBuffer;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class VariableLengthByteArrayHuffmanCodec implements BitCodec<byte[]> {
	private HuffmanCodec<Byte> byteCodec;
	private HuffmanCodec<Byte> lenCodec;
	private ByteBuffer buf;

	public VariableLengthByteArrayHuffmanCodec(byte[] alphabet, int[] freqs, byte[] lenAlphabet, int[] lenFreqs) {

		int maxLen = 0;
		for (byte len : lenAlphabet)
			if (maxLen < (0xFF & len))
				maxLen = (0xFF & len);
		
		buf = ByteBuffer.allocate(maxLen);
		HuffmanTree<Byte> tree = HuffmanCode.buildTree(freqs, Utils.autobox(alphabet));
		byteCodec = new HuffmanCodec<Byte>(tree);
		
		HuffmanTree<Byte> lenTree = HuffmanCode.buildTree(lenFreqs, Utils.autobox(lenAlphabet));
		lenCodec = new HuffmanCodec<Byte>(lenTree);
	}

	@Override
	public byte[] read(BitInputStream bis) throws IOException {
		buf.clear();
		
		Byte len = lenCodec.read (bis) ;
		
		for (int i = 0; i < (0xFF & len); i++)
			buf.put(byteCodec.read(bis));

		byte[] sequence = new byte[buf.position()];
		buf.flip();
		buf.get(sequence);
		return sequence;
	}

	@Override
	public long write(BitOutputStream bos, byte[] bytes) throws IOException {
		long len = 0;
		
		len += lenCodec.write(bos, (byte)bytes.length) ;
		
		for (byte b : bytes)
			len += byteCodec.write(bos, b);

		return len;
	}

	@Override
	public long numberOfBits(byte[] bytes) {
		long len = 0;
		
		len += lenCodec.numberOfBits((byte)bytes.length) ;
		
		for (byte b : bytes)
			len += byteCodec.numberOfBits(b);

		return len;
	}

}
