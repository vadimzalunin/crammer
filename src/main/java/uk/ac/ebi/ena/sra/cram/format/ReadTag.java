package uk.ac.ebi.ena.sra.cram.format;

public class ReadTag implements Comparable<ReadTag> {
	private static final long MAX_INT = Integer.MAX_VALUE;
	private static final long MAX_UINT = MAX_INT * 2 + 1;
	private static final long MAX_SHORT = Short.MAX_VALUE;
	private static final long MAX_USHORT = MAX_SHORT * 2 + 1;
	private static final long MAX_BYTE = Byte.MAX_VALUE;
	private static final long MAX_UBYTE = MAX_BYTE * 2 + 1;

	// non-null
	private String key;
	private String keyAndType;
	private char type;
	private Object value;

	public ReadTag(String key, char type, Object value) {
		if (key == null)
			throw new NullPointerException("Tag key cannot be null.");
		if (value == null)
			throw new NullPointerException("Tag value cannot be null.");

		this.value = value;

		if (key.length() == 2) {
			this.key = key;
			this.type = getTagValueType(value);
			keyAndType = key + ":" + getType();
		} else if (key.length() == 4) {
			this.key = key.substring(0, 2);
			this.type = key.charAt(3);
		}
	}

	public static ReadTag deriveTypeFromKeyAndType(String keyAndType, Object value) {
		if (keyAndType.length() != 4)
			throw new RuntimeException("Tag key and type must be 4 char long: " + keyAndType);

		return new ReadTag(keyAndType.substring(0, 2), keyAndType.charAt(3), value);
	}

	public static ReadTag deriveTypeFromValue(String key, Object value) {
		if (key.length() != 2)
			throw new RuntimeException("Tag key must be 2 char long: " + key);

		return new ReadTag(key, getTagValueType(value), value);
	}

	public String getKey() {
		return key;
	}

	@Override
	public int compareTo(ReadTag o) {
		return key.compareTo(o.key);
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof ReadTag))
			return false;

		ReadTag foe = (ReadTag) obj;
		if (!key.equals(foe.key))
			return false;
		if (value == null && foe.value == null)
			return true;
		if (value != null && value.equals(foe.value))
			return true;

		return false;
	}

	@Override
	public int hashCode() {
		return key.hashCode();
	}

	public Object getValue() {
		return value;
	}

	public char getType() {
		return type;
	}
	
	public String getKeyAndType() {
		return keyAndType;
	}
	
	public byte[] getValueAsByteArray () {
		return value.toString().getBytes() ;
	}
	
	public static Object restoreValueFromByteArray (char type, byte[] array) {
		switch (type) {
		case 'Z':
			return new String (array) ;
		case 'A':
			return new String (array).charAt(0) ;
		case 'f':
			return Float.valueOf(new String (array)) ;
			
		case 'I':
			return Long.valueOf(new String (array)) ;
		case 'i':
		case 'S':
			return Integer.valueOf(new String (array)) ;
		case 's':
			return Short.valueOf(new String (array)) ;
		case 'C':
			return Integer.valueOf(new String (array)) ;
		case 'c':
			return Byte.valueOf(new String (array)) ;

		default:
			throw new RuntimeException("Unknown tag type: " + type) ;
		}
	}
	
	// copied from net.sf.samtools.BinaryTagCodec 1.62:
	public static char getTagValueType(final Object value) {
		if (value instanceof String) {
			return 'Z';
		} else if (value instanceof Character) {
			return 'A';
		} else if (value instanceof Float) {
			return 'f';
		} else if (value instanceof Number) {
			if (!(value instanceof Byte || value instanceof Short || value instanceof Integer || value instanceof Long)) {
				throw new IllegalArgumentException("Unrecognized tag type " + value.getClass().getName());
			}
			return getIntegerType(((Number) value).longValue());
		} /*
		 * Note that H tag type is never written anymore, because B style is
		 * more compact. else if (value instanceof byte[]) { return 'H'; }
		 */
		else if (value instanceof byte[] || value instanceof short[] || value instanceof int[]
				|| value instanceof float[]) {
			return 'B';
		} else {
			throw new IllegalArgumentException("When writing BAM, unrecognized tag type " + value.getClass().getName());
		}
	}

	// copied from net.sf.samtools.BinaryTagCodec:
	static private char getIntegerType(final long val) {
		if (val > MAX_UINT) {
			throw new IllegalArgumentException("Integer attribute value too large to be encoded in BAM");
		}
		if (val > MAX_INT) {
			return 'I';
		}
		if (val > MAX_USHORT) {
			return 'i';
		}
		if (val > MAX_SHORT) {
			return 'S';
		}
		if (val > MAX_UBYTE) {
			return 's';
		}
		if (val > MAX_BYTE) {
			return 'C';
		}
		if (val >= Byte.MIN_VALUE) {
			return 'c';
		}
		if (val >= Short.MIN_VALUE) {
			return 's';
		}
		if (val >= Integer.MIN_VALUE) {
			return 'i';
		}
		throw new IllegalArgumentException("Integer attribute value too negative to be encoded in BAM");
	}

}
