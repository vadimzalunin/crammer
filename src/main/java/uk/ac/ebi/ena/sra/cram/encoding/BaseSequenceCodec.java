package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class BaseSequenceCodec implements BitCodec<byte[]> {
	public enum BaseCodecType {
		FLAT, RAISED, STEEP;
	}

	private BitCodec<byte[]> delegate;

	public BaseSequenceCodec(BaseCodecType type, byte[] order) {
		switch (type) {
		case FLAT:
			delegate = new BaseSequenceFlatCodec(order);
			break;
		case RAISED:
			delegate = new BaseSequenceRaisedCodec(order);
			break;
		case STEEP:
			delegate = new BaseSequenceSteepCodec(order);
			break;

		default:
			throw new IllegalArgumentException("Unknown base codec type: "
					+ type);
		}
	}

	@Override
	public byte[] read(BitInputStream bis) throws IOException {
		return delegate.read(bis);
	}

	@Override
	public long write(BitOutputStream bis, byte[] object) throws IOException {
		return delegate.write(bis, object);
	}

}
