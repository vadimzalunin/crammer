package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;
import java.util.Arrays;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

/**
 * 
 * 
 * <pre>
 * flat:
 * 1	000		0
 * 2	001		1
 * 3	010		2
 * 4	011		3
 * 5	100		4
 * 
 * normal:
 * 1	00		0
 * 2	01		1
 * 3	10		2
 * 4	110		6
 * 5	111		7
 * 
 * steep:
 * 1	0		0
 * 2	10		2
 * 3	110		6
 * 4	1110	14
 * 5	1111	15
 * </pre>
 * 
 * @author vadim
 * 
 */
class BaseRaisedCodec implements BitCodec<Byte> {

	private byte[] order;
	private int[] base2indexArray;

	public BaseRaisedCodec() {
		this("ACGTN".getBytes());
	}

	public BaseRaisedCodec(byte[] order) {
		if (order.length != 5)
			throw new IllegalArgumentException(
					"Expecting 5 bases order only but got: "
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
	public Byte read(BitInputStream bis) throws IOException {
		int bitsOneTwo = bis.readBits(2);
		if (bitsOneTwo < 3)
			return order[bitsOneTwo];

		if (bis.readBit())
			return order[4];

		return order[3];
	}

	@Override
	public long write(BitOutputStream bis, Byte base) throws IOException {
		int index = base2indexArray[base];
		if (index < 0 || index > order.length)
			throw new IllegalArgumentException("Invalid base byte: " + base);

		switch (index) {
		case 0:
			bis.write(0, 2);
			return 2;
		case 1:
			bis.write(1, 2);
			return 2;
		case 2:
			bis.write(2, 2);
			return 2;
		case 3:
			bis.write(6, 3);
			return 3;
		case 4:
			bis.write(7, 3);
			return 3;

		default:
			throw new RuntimeException("Unkown base: " + base);
		}

	}

	@Override
	public long numberOfBits(Byte base) {
		int index = base2indexArray[base];
		if (index < 0 || index > order.length)
			throw new IllegalArgumentException("Invalid base byte: " + base);

		if (index < 3)
			return 2L;
		else
			return 3L;
	}
}
