package uk.ac.ebi.ena.sra.cram;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.format.CramFormatException;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

public class CramTools {
	public static final String BAM_STATS_COMMAND = "stats";
	public static final String CRAM2BAM_COMMAND = "bam";
	public static final String BAM2CRAM_COMMAND = "cram";
	public static final String CRAM2TEXT_COMMAND = "view";

	public static void main(String[] args) throws IOException,
			CramFormatException {
		Params params = new Params();
		JCommander jc = new JCommander(params);
		jc.setProgramName("cramtools");

		BamStats.Params bamStatsParams = new BamStats.Params();
		Cram2Bam.Params cram2BamParams = new Cram2Bam.Params();
		Bam2Cram.Params bam2CramParams = new Bam2Cram.Params();
		Cram2Text.Params cram2TextParams = new Cram2Text.Params();

		jc.addCommand(BAM_STATS_COMMAND, bamStatsParams);
		jc.addCommand(CRAM2BAM_COMMAND, cram2BamParams);
		jc.addCommand(BAM2CRAM_COMMAND, bam2CramParams);
		jc.addCommand(CRAM2TEXT_COMMAND, cram2TextParams);

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

		if (BAM_STATS_COMMAND.equals(command))
			BamStats.main(commandArgs);
		else if (CRAM2BAM_COMMAND.equals(command))
			Cram2Bam.main(commandArgs);
		if (BAM2CRAM_COMMAND.equals(command))
			Bam2Cram.main(commandArgs);
		if (CRAM2TEXT_COMMAND.equals(command))
			Cram2Text.main(commandArgs);

	}

	private static class Params {
		@Parameter(names = { "-h", "--help" }, description = "Print help and quit")
		private boolean help = false;
	}
}
