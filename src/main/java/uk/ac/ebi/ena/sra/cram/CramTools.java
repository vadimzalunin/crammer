package uk.ac.ebi.ena.sra.cram;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

public class CramTools {
	public static final String CRAM2BAM_COMMAND = "bam";
	public static final String BAM2CRAM_COMMAND = "cram";
	public static final String CRAMINDEX_COMMAND = "index";

	public static void main(String[] args) throws Exception {

		Params params = new Params();
		JCommander jc = new JCommander(params);
		jc.setProgramName("cramtools");

		Cram2Bam.Params cram2BamParams = new Cram2Bam.Params();
		Bam2Cram.Params bam2CramParams = new Bam2Cram.Params();
		CramIndexer.Params cramIndexerParams = new CramIndexer.Params();

		jc.addCommand(CRAM2BAM_COMMAND, cram2BamParams);
		jc.addCommand(BAM2CRAM_COMMAND, bam2CramParams);
		jc.addCommand(CRAMINDEX_COMMAND, cramIndexerParams);

		jc.parse(args);

		if (params.logLevel != null) {
			Logger.getRootLogger().setLevel(params.logLevel);
			String[] newArgs = new String[args.length - 2];
			System.arraycopy(args, 2, newArgs, 0, newArgs.length);
			args = newArgs;
		}

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
			Bam2Cram.main(commandArgs);
		else if (CRAMINDEX_COMMAND.equals(command))
			CramIndexer.main(commandArgs);

	}

	@Parameters(commandDescription = "CRAM tools. Version 0.65")
	private static class Params  {
		@Parameter(names = { "-h", "--help" }, description = "Print help and quit")
		private boolean help = false;

		@Parameter(names = { "-l", "--log-level" }, description = "Change log level", converter = LevelConverter.class)
		private Level logLevel = null;
	}

	public static class LevelConverter implements IStringConverter<Level> {

		@Override
		public Level convert(String s) {
			return Level.toLevel(s);
		}

	}
}
