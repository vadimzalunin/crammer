package net.sf.cram.structure;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import net.sf.cram.EncodingID;
import net.sf.cram.EncodingKey;
import net.sf.cram.EncodingParams;
import net.sf.cram.encoding.NullEncoding;
import net.sf.cram.io.ByteBufferUtils;
import net.sf.picard.util.Log;

public class CompressionHeader {
	private static Log log = Log.getInstance(CompressionHeader.class);

	public boolean readNamesIncluded;
	public boolean AP_seriesDelta=true;

	public Map<EncodingKey, EncodingParams> eMap;
	public Map<Integer, EncodingParams> tMap;

	public SubstitutionMatrix substitutionMatrix ;

	public List<Integer> externalIds;

	public byte[][][] dictionary;

	public CompressionHeader() {
	}

	public CompressionHeader(InputStream is) throws IOException {
		read(is);
	}

	private byte[][][] parseDictionary(byte[] bytes) {
		List<List<byte[]>> dictionary = new ArrayList<List<byte[]>>();
		{
			int i = 0;
			while (i < bytes.length) {
				List<byte[]> list = new ArrayList<byte[]>();
				while (bytes[i] != 0) {
					list.add(Arrays.copyOfRange(bytes, i, i + 3));
					i += 3;
				}
				i++;
				dictionary.add(list);
			}
		}

		int maxWidth = 0;
		for (List<byte[]> list : dictionary)
			maxWidth = Math.max(maxWidth, list.size());

		byte[][][] array = new byte[dictionary.size()][][];
		for (int i = 0; i < dictionary.size(); i++) {
			List<byte[]> list = dictionary.get(i);
			array[i] = (byte[][]) list.toArray(new byte[list.size()][]);
		}

		return array;
	}

	private byte[] dictionaryToByteArray() {
		int size = 0;
		for (int i = 0; i < dictionary.length; i++) {
			for (int j = 0; j < dictionary[i].length; j++)
				size += dictionary[i][j].length;
			size++;
		}

		byte[] bytes = new byte[size];
		ByteBuffer buf = ByteBuffer.wrap(bytes);
		for (int i = 0; i < dictionary.length; i++) {
			for (int j = 0; j < dictionary[i].length; j++)
				buf.put(dictionary[i][j]);
			buf.put((byte) 0);
		}

		return bytes;
	}

	public byte[][] getTagIds(int id) {
		return dictionary[id];
	}

	public void read(byte[] data) {
		ByteArrayInputStream bais = new ByteArrayInputStream(data);
		try {
			read(bais);
		} catch (IOException e) {
			throw new RuntimeException("This should have never happened.");
		}
	}

	public void read(InputStream is) throws IOException {
		{ // preservation map:
			int byteSize = ByteBufferUtils.readUnsignedITF8(is);
			byte[] bytes = new byte[byteSize];
			ByteBufferUtils.readFully(bytes, is);
			ByteBuffer buf = ByteBuffer.wrap(bytes);

			int mapSize = ByteBufferUtils.readUnsignedITF8(buf);
			for (int i = 0; i < mapSize; i++) {
				String key = new String(new byte[] { buf.get(), buf.get() });
				if ("RN".equals(key))
					readNamesIncluded = buf.get() == 1 ? true : false;
				else if ("AP".equals(key))
					AP_seriesDelta = buf.get() == 1 ? true : false;
				else if ("TD".equals(key)) {
					int size = ByteBufferUtils.readUnsignedITF8(buf);
					byte[] dictionaryBytes = new byte[size];
					buf.get(dictionaryBytes);
					dictionary = parseDictionary(dictionaryBytes);
				} else if ("SM".equals(key)) {
					// parse subs matrix here:
					byte[] matrixBytes = new byte[5] ;
					buf.get(matrixBytes);
					substitutionMatrix = new SubstitutionMatrix(matrixBytes) ;
				} else
					throw new RuntimeException("Unknown preservation map key: "
							+ key);
			}
		}

		{ // encoding map:
			int byteSize = ByteBufferUtils.readUnsignedITF8(is);
			byte[] bytes = new byte[byteSize];
			ByteBufferUtils.readFully(bytes, is);
			ByteBuffer buf = ByteBuffer.wrap(bytes);

			int mapSize = ByteBufferUtils.readUnsignedITF8(buf);
			eMap = new TreeMap<EncodingKey, EncodingParams>();
			for (EncodingKey key : EncodingKey.values())
				eMap.put(key, NullEncoding.toParam());

			for (int i = 0; i < mapSize; i++) {
				String key = new String(new byte[] { buf.get(), buf.get() });
				EncodingKey eKey = EncodingKey.byFirstTwoChars(key);
				if (eKey == null)
					throw new RuntimeException("Unknown encoding key: " + key);

				EncodingID id = EncodingID.values()[buf.get()];
				int paramLen = ByteBufferUtils.readUnsignedITF8(buf);
				byte[] paramBytes = new byte[paramLen];
				buf.get(paramBytes);

				eMap.put(eKey, new EncodingParams(id, paramBytes));

				log.debug(String.format("FOUND ENCODING: %s, %s, %s.",
						eKey.name(), id.name(),
						Arrays.toString(Arrays.copyOf(paramBytes, 20))));
			}
		}

		{ // tag encoding map:
			int byteSize = ByteBufferUtils.readUnsignedITF8(is);
			byte[] bytes = new byte[byteSize];
			ByteBufferUtils.readFully(bytes, is);
			ByteBuffer buf = ByteBuffer.wrap(bytes);

			int mapSize = ByteBufferUtils.readUnsignedITF8(buf);
			tMap = new TreeMap<Integer, EncodingParams>();
			for (int i = 0; i < mapSize; i++) {
				int key = ByteBufferUtils.readUnsignedITF8(buf);

				EncodingID id = EncodingID.values()[buf.get()];
				int paramLen = ByteBufferUtils.readUnsignedITF8(buf);
				byte[] paramBytes = new byte[paramLen];
				buf.get(paramBytes);

				tMap.put(key, new EncodingParams(id, paramBytes));
			}
		}
	}

	public byte[] toByteArray() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		write(baos);
		return baos.toByteArray();
	}

	public void write(OutputStream os) throws IOException {

		{ // preservation map:
			ByteBuffer mapBuf = ByteBuffer.allocate(1024 * 100);
			ByteBufferUtils.writeUnsignedITF8(4, mapBuf);

			mapBuf.put("RN".getBytes());
			mapBuf.put((byte) (readNamesIncluded ? 1 : 0));

			mapBuf.put("AP".getBytes());
			mapBuf.put((byte) (AP_seriesDelta ? 1 : 0));

			mapBuf.put("SM".getBytes());
			mapBuf.put(substitutionMatrix.getEncodedMatrix());

			mapBuf.put("TD".getBytes());
			{

				byte[] dBytes = dictionaryToByteArray();
				ByteBufferUtils.writeUnsignedITF8(dBytes.length, mapBuf);
				mapBuf.put(dBytes);
			}

			mapBuf.flip();
			byte[] mapBytes = new byte[mapBuf.limit()];
			mapBuf.get(mapBytes);

			ByteBufferUtils.writeUnsignedITF8(mapBytes.length, os);
			os.write(mapBytes);
		}

		{ // encoding map:
			int size = 0;
			for (EncodingKey eKey : eMap.keySet()) {
				if (eMap.get(eKey).id != EncodingID.NULL)
					size++;
			}

			ByteBuffer mapBuf = ByteBuffer.allocate(1024 * 100);
			ByteBufferUtils.writeUnsignedITF8(size, mapBuf);
			for (EncodingKey eKey : eMap.keySet()) {
				if (eMap.get(eKey).id == EncodingID.NULL)
					continue;

				mapBuf.put((byte) eKey.name().charAt(0));
				mapBuf.put((byte) eKey.name().charAt(1));

				EncodingParams params = eMap.get(eKey);
				mapBuf.put((byte) (0xFF & params.id.ordinal()));
				ByteBufferUtils.writeUnsignedITF8(params.params.length, mapBuf);
				mapBuf.put(params.params);
			}
			mapBuf.flip();
			byte[] mapBytes = new byte[mapBuf.limit()];
			mapBuf.get(mapBytes);

			ByteBufferUtils.writeUnsignedITF8(mapBytes.length, os);
			os.write(mapBytes);
		}

		{ // tag encoding map:
			ByteBuffer mapBuf = ByteBuffer.allocate(1024 * 100);
			ByteBufferUtils.writeUnsignedITF8(tMap.size(), mapBuf);
			for (Integer eKey : tMap.keySet()) {
				ByteBufferUtils.writeUnsignedITF8(eKey, mapBuf);

				EncodingParams params = tMap.get(eKey);
				mapBuf.put((byte) (0xFF & params.id.ordinal()));
				ByteBufferUtils.writeUnsignedITF8(params.params.length, mapBuf);
				mapBuf.put(params.params);
			}
			mapBuf.flip();
			byte[] mapBytes = new byte[mapBuf.limit()];
			mapBuf.get(mapBytes);

			ByteBufferUtils.writeUnsignedITF8(mapBytes.length, os);
			os.write(mapBytes);
		}
	}

}
