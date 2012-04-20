package uk.ac.ebi.ena.sra.cram;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class CramBatch {

	public static void main(String[] args) throws Exception {
		List<String> sources = new ArrayList<String>();
		if (args.length == 0)
			sources.add("cram.options");
		else
			sources.add(args[0]);

		for (String sourceFileName : sources) {
			File optionsFile = new File(sourceFileName);
			InputStream is;
			if (!optionsFile.exists())
				is = System.in;
			else
				is = new FileInputStream(optionsFile);
			Scanner scanner = new Scanner(is);
			String[] options = scanner.nextLine().split(" ");
			CramTools.main(options);
		}
	}
}
