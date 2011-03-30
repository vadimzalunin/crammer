package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.Arrays;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

class BaseSequenceFlatCodec implements BitCodec<byte[]> {
	private byte[] order = "ACGTNS".getBytes();
	private int[] base2indexArray;
	private ByteBuffer buffer = ByteBuffer.allocate(1024);

	public BaseSequenceFlatCodec(byte[] order) {
		if (order.length != 6)
			throw new IllegalArgumentException(
					"Expecting 5 bases and 1 stop only but got: "
							+ new String(order));

		this.order = order;
		this.base2indexArray = new int[255];
		Arrays.fill(base2indexArray, -1);
		for (int i = 0; i < 255; i++)
			base2indexArray[i] = -1;

		int index = 0;
		for (byte base : order)
			base2indexArray[base] = index++;
	}

	@Override
	public byte[] read(BitInputStream bis) throws IOException {
		buffer.clear();

		while (true) {
			int threeBits = bis.readBits(3);

			if (threeBits >= order.length)
				throw new RuntimeException("Unexpected base flat code: "
						+ threeBits);

			if (order[threeBits] == 'S')
				break;
			buffer.put(order[threeBits]);
		}

		buffer.flip() ;
		byte[] seq = new byte[buffer.limit()];
		buffer.get(seq);
		return seq;
	}

	@Override
	public long write(BitOutputStream bis, byte[] bases) throws IOException {
		for (byte base : bases)
			bis.writeBits(base2indexArray[base], 3);

		bis.writeBits(base2indexArray['S'], 3);
		return (bases.length + 1) * 3;
	}
}
