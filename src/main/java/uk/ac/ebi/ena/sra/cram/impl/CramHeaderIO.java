package uk.ac.ebi.ena.sra.cram.impl;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;

import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;

public class CramHeaderIO {
	private static final String MAGIC = "CRAM";

	public static void write(CramHeader header, OutputStream os)
			throws IOException {
		DataOutputStream headerDOS = new DataOutputStream(os);
		headerDOS.writeUTF(MAGIC);
		headerDOS.writeUTF(header.getVersion());

		headerDOS.writeInt(header.getReferenceSequences().size());
		for (CramReferenceSequence s : header.getReferenceSequences()) {
			headerDOS.writeUTF(s.getName());
			headerDOS.writeInt(s.getLength());
		}
		headerDOS.flush();
	}

	public static CramHeader read(InputStream is) throws IOException {
		CramHeader header = new CramHeader();
		DataInputStream dis = new DataInputStream(is);
		String magick = dis.readUTF();
		if (!MAGIC.equals(magick))
			throw new RuntimeException("Not recognized as CRAM format.");
		header.setVersion(dis.readUTF());
		int seqCount = dis.readInt();
		ArrayList<CramReferenceSequence> list = new ArrayList<CramReferenceSequence>();
		for (int i = 0; i < seqCount; i++) {
			CramReferenceSequence cramSeq = new CramReferenceSequence();
			cramSeq.setName(dis.readUTF());
			cramSeq.setLength(dis.readInt());
			list.add(cramSeq);
		}
		header.setReferenceSequences(list);
		return header;
	}
}
