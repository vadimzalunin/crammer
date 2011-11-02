package net.sf.samtools;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.util.SeekableStream;

public class ReferenceDiscovery {
	public static Map<Object, ReferenceSequenceFile> referenceFactory = new HashMap<Object, ReferenceSequenceFile>();

	public static ReferenceSequenceFile probeLocation(String location) {
		ReferenceSequenceFile referenceSequenceFile = referenceFactory.get(location);
		if (referenceSequenceFile != null)
			return referenceSequenceFile;

		String baseName = location.replaceFirst(".cram$", "");

		File file = new File(baseName + ".fa");
		if (file.exists())
			return ReferenceSequenceFileFactory.getReferenceSequenceFile(file);

		file = new File(baseName + ".ref");
		if (file.exists()) {
			Scanner scanner = null;
			try {
				scanner = new Scanner(file);
			} catch (FileNotFoundException e) {
				throw new RuntimeException(e);
			}
			if (scanner.hasNextLine()) {
				file = new File(scanner.nextLine());
				if (file.exists()) {
					return ReferenceSequenceFileFactory.getReferenceSequenceFile(file);
				}
			}
		}
		return null;
	}

	public static ReferenceSequenceFile findReferenceSequenceFileOrFail(Object... sources) {
		ReferenceSequenceFile referenceSequenceFile = null;

		if (sources != null)
			for (Object source : sources) {
				if (source == null)
					continue;
				ReferenceSequenceFile foundFile = referenceFactory.get(source);
				if (foundFile != null)
					return foundFile;

				// last resort, using a blind-walk method in a black room to
				// find a
				// Schrödinger's cat
				// which is probably not there:
				String name = null;
				if (source instanceof String)
					name = (String) source;
				else if (source instanceof File)
					name = ((File) source).getAbsolutePath();
				else if (source instanceof SeekableStream)
					name = ((SeekableStream) source).getSource();

				if (name == null)
					continue ;

				referenceSequenceFile = probeLocation(name);
				if (referenceSequenceFile != null)
					return referenceSequenceFile;

				// could be an index file:
				referenceSequenceFile = probeLocation(name.replaceAll(".crai$", ""));
				if (referenceSequenceFile != null)
					return referenceSequenceFile;
			}

		// try default rsf from java properties:
		String refProperty = System.getProperty("reference");
		if (refProperty != null) {
			File file = new File(refProperty);
			if (file.isFile())
				return ReferenceSequenceFileFactory.getReferenceSequenceFile(file);
		}

		throw new RuntimeException("Reference sequence file not found.");
	}
}
