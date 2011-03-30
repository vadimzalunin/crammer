package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class BaseCodec implements BitCodec<Byte> {
	public enum BaseCodecType {
		FLAT, RAISED, STEEP;
	}

	private BitCodec<Byte> delegate;

	public BaseCodec(BaseCodecType type, byte[] order) {
		switch (type) {
		case FLAT:
			delegate = new BaseFlatCodec(order);
			break;
		case RAISED:
			delegate = new BaseRaisedCodec(order);
			break;
		case STEEP:
			delegate = new BaseSteepCodec(order);
			break;

		default:
			throw new IllegalArgumentException("Unknown base codec type: "
					+ type);
		}
	}

	@Override
	public Byte read(BitInputStream bis) throws IOException {
		return delegate.read(bis);
	}

	@Override
	public long write(BitOutputStream bis, Byte object) throws IOException {
		return delegate.write(bis, object);
	}

}
