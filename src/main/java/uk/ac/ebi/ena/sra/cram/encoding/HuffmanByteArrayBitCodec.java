package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.ByteArrayOutputStream;
import java.io.IOException;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class HuffmanByteArrayBitCodec implements ByteArrayBitCodec {
	private HuffmanByteCodec2 byteCodec;
	private ByteArrayOutputStream baos = new ByteArrayOutputStream(1024);

	private final Stats stats = new Stats();
	private String name;

	public HuffmanByteArrayBitCodec(byte[] alphabet, int[] freqs) {
		HuffmanTree<Byte> tree = HuffmanCode.buildTree(freqs, Utils.autobox(alphabet));
		byteCodec = new HuffmanByteCodec2(tree);
	}

	@Override
	public byte[] read(BitInputStream bis, int len) throws IOException {
		stats.arraysRead++;
		stats.bytesRead += len;
		baos.reset();

		for (int i = 0; i < len; i++)
			baos.write(byteCodec.read(bis));

		return baos.toByteArray();
	}

	@Override
	public long write(BitOutputStream bos, byte[] bytes) throws IOException {
		stats.arraysWritten++ ;
		stats.bytesWritten += bytes.length ;

		long len = 0;
		for (byte b : bytes)
			len += byteCodec.write(bos, b);
		
		stats.nofBis += len ;

		return len;
	}

	@Override
	public long numberOfBits(byte[] bytes) {
		long bitCount = 0;
		for (byte b : bytes)
			bitCount += byteCodec.numberOfBits(b);

		return bitCount;
	}

	@Override
	public Stats getStats() {
		return stats;
	}

	@Override
	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}
	
	@Override
	public String toString() {
		return String
				.format("%s: written objects=%d; written bits=%d; written bits per object=%.2f",
						name, stats.bytesWritten, stats.nofBis, (double) stats.nofBis/stats.bytesWritten);
	}
	
}
