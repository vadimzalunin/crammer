package uk.ac.ebi.ena.sra.cram;

import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

public class CramTools {
	public static final String CRAM2BAM_COMMAND = "bam";
	public static final String BAM2CRAM_COMMAND = "cram";

	private static void setupLogger() {
		PatternLayout layout = new PatternLayout(
				"%d{ABSOLUTE} %5p %c{1}:%L - %m%n");

		ConsoleAppender appender = new ConsoleAppender(layout, "System.err");
		appender.setThreshold(Level.INFO);

		Logger.getRootLogger().addAppender(appender);
		Logger.getRootLogger().setLevel(Level.ALL);
	}

	public static void main(String[] args) throws Exception {
		setupLogger();

		Params params = new Params();
		JCommander jc = new JCommander(params);
		jc.setProgramName("cramtools");

		BamStats.Params bamStatsParams = new BamStats.Params();
		Cram2Bam.Params cram2BamParams = new Cram2Bam.Params();
		CramWriterReadTest.Params bam2CramParams = new CramWriterReadTest.Params();

		jc.addCommand(CRAM2BAM_COMMAND, cram2BamParams);
		jc.addCommand(BAM2CRAM_COMMAND, bam2CramParams);

		jc.parse(args);
		String command = jc.getParsedCommand();

		if (command == null || params.help) {
			StringBuilder sb = new StringBuilder();
			sb.append("\n");
			jc.usage(sb);

			System.out.println(sb.toString());
			return;
		}

		String[] commandArgs = new String[args.length - 1];
		System.arraycopy(args, 1, commandArgs, 0, commandArgs.length);

		if (CRAM2BAM_COMMAND.equals(command))
			Cram2Bam.main(commandArgs);
		else if (BAM2CRAM_COMMAND.equals(command))
			CramWriterReadTest.main(commandArgs);

	}

	private static class Params {
		@Parameter(names = { "-h", "--help" }, description = "Print help and quit")
		private boolean help = false;
	}
}
