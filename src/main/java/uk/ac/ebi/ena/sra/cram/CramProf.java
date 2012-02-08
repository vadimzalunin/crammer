package uk.ac.ebi.ena.sra.cram;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.Scanner;

public class CramProf {

	public static void main(String[] args) throws Exception {
		File optionsFile = new File("options");
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
