package net.sf.cram.encoding.read_features;

class BaseTransitions {
	public static byte getBaseForTransition(byte reference, int transition) {
		switch (reference) {
		case 'A':
		case 'a':
			switch (transition) {
			case 0:
				return 'C';
			case 1:
				return 'G';
			case 2:
				return 'T';
			case 3:
				return 'N';

			default:
				throw new IllegalArgumentException(
						String.format(
								"Unknown transition when restoring base from transition: %c, %d",
								(char) reference, transition));
			}

		case 'C':
		case 'c':
			switch (transition) {
			case 0:
				return 'A';
			case 1:
				return 'G';
			case 2:
				return 'T';
			case 3:
				return 'N';

			default:
				throw new IllegalArgumentException(
						String.format(
								"Unknown transition when restoring base from transition: %c, %d",
								(char) reference, transition));
			}
		case 'G':
		case 'g':
			switch (transition) {
			case 0:
				return 'A';
			case 1:
				return 'C';
			case 2:
				return 'T';
			case 3:
				return 'N';

			default:
				throw new IllegalArgumentException(
						String.format(
								"Unknown transition when restoring base from transition: %c, %d",
								(char) reference, transition));
			}
		case 'T':
		case 't':
			switch (transition) {
			case 0:
				return 'A';
			case 1:
				return 'C';
			case 2:
				return 'G';
			case 3:
				return 'N';

			default:
				throw new IllegalArgumentException(
						String.format(
								"Unknown transition when restoring base from transition: %c, %d",
								(char) reference, transition));
			}
			// case 'N':
		default:
			switch (transition) {
			case 0:
				return 'A';
			case 1:
				return 'C';
			case 2:
				return 'G';
			case 3:
				return 'T';

			default:
				throw new IllegalArgumentException(
						String.format(
								"Unknown transition when restoring base from transition: %c, %d",
								(char) reference, transition));
			}

			// default:
			// throw new IllegalArgumentException(
			// String.format(
			// "Unknown reference base when restoring base from transition: %c, %d",
			// (char) reference, transition));
		}
	}

	public static int getBaseTransition(byte from, byte to) {
		switch (from) {
		case 'A':
		case 'a':
			switch (to) {
			case 'C':
				return 0;
			case 'G':
				return 1;
			case 'T':
				return 2;
			case 'N':
				return 3;

			default:
				throw new IllegalArgumentException(
						"Unknown base when calcualting base transition: "
								+ (char) from);
			}

		case 'C':
		case 'c':
			switch (to) {
			case 'A':
				return 0;
			case 'G':
				return 1;
			case 'T':
				return 2;
			case 'N':
				return 3;

			default:
				throw new IllegalArgumentException(
						"Unknown base when calcualting base transition: "
								+ (char) from);
			}
		case 'G':
		case 'g':
			switch (to) {
			case 'A':
				return 0;
			case 'C':
				return 1;
			case 'T':
				return 2;
			case 'N':
				return 3;

			default:
				throw new IllegalArgumentException(
						"Unknown base when calcualting base transition: "
								+ (char) from);
			}
		case 'T':
		case 't':
			switch (to) {
			case 'A':
				return 0;
			case 'C':
				return 1;
			case 'G':
				return 2;
			case 'N':
				return 3;

			default:
				throw new IllegalArgumentException(
						"Unknown base when calcualting base transition: "
								+ (char) from);
			}
			// case 'N':
		default:
			switch (to) {
			case 'A':
				return 0;
			case 'C':
				return 1;
			case 'G':
				return 2;
			case 'T':
				return 3;

			default:
				throw new IllegalArgumentException(
						"Unknown base when calcualting base transition: "
								+ (char) from);
			}

			// default:
			// throw new IllegalArgumentException(
			// "Unknown base when calcualting base transition: "
			// + (char) from);
		}
	}
}
