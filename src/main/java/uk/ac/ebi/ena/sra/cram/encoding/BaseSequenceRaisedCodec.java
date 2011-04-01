package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;
import java.nio.ByteBuffer;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

class BaseSequenceRaisedCodec implements BitCodec<byte[]> {
	private byte[] order = "ACGTNS".getBytes();
	private BitCode[] codes;
	private ByteBuffer buffer = ByteBuffer.allocate(1024);

	public BaseSequenceRaisedCodec(byte[] order) {
		if (order.length != 6)
			throw new IllegalArgumentException(
					"Expecting 5 bases and 1 stop only but got: "
							+ new String(order));

		this.order = order;
		this.codes = new BitCode[255];

		codes[order[0]] = new BitCode(0, 2);
		codes[order[1]] = new BitCode(1, 2);
		codes[order[2]] = new BitCode(2, 2);
		codes[order[3]] = new BitCode(6, 3);
		codes[order[4]] = new BitCode(14, 4);
		codes[order[5]] = new BitCode(15, 4);
	}

	@Override
	public byte[] read(BitInputStream bis) throws IOException {
		buffer.clear();

		while (true) {
			int bitsOneTwo = bis.readBits(2);
			byte base;
			switch (bitsOneTwo) {
			case 0:
			case 1:
			case 2:
				base = order[bitsOneTwo];
				break;

			default:
				if (!bis.readBit())
					base = order[3];
				else {
					if (!bis.readBit())
						base = order[4];
					else
						base = order[5];
				}
				break;
			}
			if (base == 'S')
				break;
			buffer.put(base);
		}

		buffer.flip();
		byte[] seq = new byte[buffer.limit()];
		buffer.get(seq);
		return seq;
	}

	@Override
	public long write(BitOutputStream bos, byte[] bases) throws IOException {
		int length = 0;
		for (byte base : bases) {
			if (base == 'S')
				throw new IllegalArgumentException(
						"Stop code is not allowed in a base sequence.");
			BitCode code = codes[base];
			bos.write(code.value, code.bits);
			length += code.bits;
		}

		bos.write(codes['S'].value, codes['S'].bits);
		length += codes['S'].bits;

		return length;
	}

	@Override
	public long numberOfBits(byte[] bases) {
		int length = 0;
		for (byte base : bases) {
			if (base == 'S')
				throw new IllegalArgumentException(
						"Stop code is not allowed in a base sequence.");
			BitCode code = codes[base];
			length += code.bits;
		}

		length += codes['S'].bits;

		return length;
	}
}
