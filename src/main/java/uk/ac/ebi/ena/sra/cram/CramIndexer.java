package uk.ac.ebi.ena.sra.cram;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.ReferenceDiscovery;
import net.sf.samtools.SAMRecord;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;
import uk.ac.ebi.ena.sra.cram.format.text.CramRecordFormat;
import uk.ac.ebi.ena.sra.cram.impl.ByteArraySequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.impl.CramHeaderIO;
import uk.ac.ebi.ena.sra.cram.impl.SequentialCramReader;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.FileConverter;

public class CramIndexer {
	private static Logger log = Logger.getLogger(CramIndexer.class);

	private static InputStream createCramInputStream(File file) throws IOException {
		FileInputStream fis = new FileInputStream(file);

		// gzip magic:
		if (fis.read() == 31 && fis.read() == 139)
			return new GZIPInputStream(new BufferedInputStream(new FileInputStream(file)));

		return new BufferedInputStream(new FileInputStream(file));
	}

	public static void main(String[] args) throws Exception {
		Params params = new Params();
		JCommander jc = new JCommander(params);
		jc.parse(args);

		if (args.length == 0 || params.help) {
			StringBuilder sb = new StringBuilder();
			sb.append("\n");
			jc.usage(sb);

			System.out.println("Version " + CramIndexer.class.getPackage().getImplementationVersion());
			System.out.println(sb.toString());
			return;
		}

		ReferenceSequenceFile referenceSequenceFile;
		if (params.reference == null)
			referenceSequenceFile = ReferenceDiscovery.findReferenceSequenceFileOrFail(params.cramFile);
		else
			referenceSequenceFile = Utils.createIndexedFastaSequenceFile(params.reference);

		InputStream cramIS = createCramInputStream(params.cramFile);

		OutputStream indexOS = null;
		if (params.indexFile == null) {
			indexOS = new GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(params.cramFile + ".crai")));
		} else
			indexOS = new GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(params.indexFile)));

		index(referenceSequenceFile, cramIS, indexOS, params.resolution);
		indexOS.close();
	}

	public static class CountingInputStream extends InputStream {
		private InputStream delegate;
		private long position = 0;
		private boolean debug = false;

		public void setDebug(boolean debug) {
			this.debug = debug;
		}

		public CountingInputStream(InputStream inputStream) {
			delegate = inputStream;
		}

		@Override
		public int read() throws IOException {
			position++;
			int read = delegate.read();
			if (debug)
				System.out.printf("pos=%d\tread=%d\n", position, read);
			return read;
		}

		public int read(byte[] b) throws IOException {
			int read = delegate.read(b);
			position += read;
			return read;
		}

		public int read(byte[] b, int off, int len) throws IOException {
			int read = delegate.read(b, off, len);
			position += read;
			return read;
		}

		public long skip(long n) throws IOException {
			long skipped = delegate.skip(n);
			position += skipped;
			return skipped;
		}

		public int available() throws IOException {
			return delegate.available();
		}

		public void close() throws IOException {
			delegate.close();
		}

		public void mark(int readlimit) {
			delegate.mark(readlimit);
		}

		public void reset() throws IOException {
			delegate.reset();
			position = 0;
		}

		public boolean markSupported() {
			return delegate.markSupported();
		}

		public long getPosition() {
			return position;
		}
	}

	public static void index(File cramFile, long resolution) throws Exception {
		ReferenceSequenceFile referenceSequenceFile;
		referenceSequenceFile = ReferenceDiscovery.findReferenceSequenceFileOrFail(cramFile);

		InputStream cramIS = new BufferedInputStream(new FileInputStream(cramFile));

		OutputStream indexOS = null;
		indexOS = new GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(cramFile + ".crai")));

		index(referenceSequenceFile, cramIS, indexOS, resolution);
		indexOS.close();
	}

	public static void index(ReferenceSequenceFile referenceSequenceFile, InputStream cramInputStream,
			OutputStream indexOS, long resolution) throws Exception {

		Utils.isCRAM(cramInputStream);

		PrintStream indexPS = new PrintStream(indexOS);

		CountingInputStream cis = new CountingInputStream(cramInputStream);
		DataInputStream cramDIS = new DataInputStream(cis);
		CramHeader cramHeader = CramHeaderIO.read(Utils.getNextChunk(cramDIS));
		log.info("CRAM format version: " + cramHeader.getVersion()) ;
		long firstBlockStart = cis.getPosition();
		int blockNumber = 0;
		List<String> sequences = new ArrayList<String>();
		for (CramReferenceSequence seq : cramHeader.getReferenceSequences()) {
			sequences.add(seq.getName());
		}

		long time1 = System.currentTimeMillis();
		CramRecordFormat cramRecordFormat = new CramRecordFormat();

		String prevSeqName = null;

		List<TreeMap<Long, Long>> indexList = new ArrayList<TreeMap<Long, Long>>();
		List<Long> alPositions = new ArrayList<Long>();
		List<Long> recordNumbersInBlock = new ArrayList<Long>();
		List<Long> byteOffsets = new ArrayList<Long>();
		List<Byte> bitOffsets = new ArrayList<Byte>();
		List<Integer> blockNumbers = new ArrayList<Integer>();
		List<Long> blockStarts = new ArrayList<Long>();
		List<Integer> sequenceNumbers = new ArrayList<Integer>();

		int blockCounter = 0;
		ByteArraySequenceBaseProvider provider = null;

		while (true) {
			long blockTime1 = System.currentTimeMillis();
			long blockStart = cis.getPosition() + 4;
			DataInputStream uncompressedDIS = null;
			CountingInputStream uncompressedCIS = null;

			uncompressedDIS = Utils.getNextChunk(cramDIS);
			if (uncompressedDIS == null)
				break;

			uncompressedCIS = new CountingInputStream(uncompressedDIS);

			SequentialCramReader reader = new SequentialCramReader(new DataInputStream(uncompressedCIS), null,
					cramHeader);
			CramRecordBlock readBlock = reader.readBlock();
			if (readBlock == null)
				break;

			int sequenceIndex = sequences.indexOf(readBlock.getSequenceName());
			TreeMap<Long, Long> map;
			if (indexList.size() < sequenceIndex + 1) {
				map = new TreeMap<Long, Long>();
				indexList.add(map);
			} else
				map = indexList.get(sequenceIndex);

			sequenceNumbers.add(sequenceIndex);

			cramRecordFormat.setSequenceID(readBlock.getSequenceName());

			String seqName = readBlock.getSequenceName();
			if (SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(seqName))
				// that's it folks!
				break;

			if (prevSeqName == null)
				prevSeqName = seqName;

			prevSeqName = seqName;

			if (provider == null || !seqName.equals(prevSeqName)) {
				log.debug("Starting sequence: " + seqName);
				byte[] refBases = Utils.getReferenceSequenceBases(referenceSequenceFile, seqName);
				provider = new ByteArraySequenceBaseProvider(refBases);
			}
			reader.setReferenceBaseProvider(provider);

			CramRecord cramRecord = null;

			long blockEnd = uncompressedCIS.getPosition();
			for (long i = 0; i < readBlock.getRecordCount(); i++) {
				cis.setDebug(false);
				try {
					if (i % resolution == 0) {
						cis.setDebug(false);
						long byteOffset = uncompressedCIS.getPosition() - blockEnd;
						byte bitOffset = reader.getHangingBits();
						cramRecord = reader.readRecord();

						sequenceNumbers.add(sequenceIndex);
						blockNumbers.add(blockNumber);
						recordNumbersInBlock.add(i);
						byteOffsets.add(byteOffset);
						bitOffsets.add(bitOffset);
						alPositions.add(cramRecord.getAlignmentStart());

						// System.out.println(cramRecord.toString());
						String line = String.format("%s\t%d\t%d\t%d\t%d\t%d\t%d\n", readBlock.getSequenceName(),
								blockNumber, blockStart, i, byteOffset, bitOffset, cramRecord.getAlignmentStart());
						indexPS.printf(line);
						// System.out.println(line);
						// System.out.printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
						// readBlock.getSequenceName(), blockNumber,
						// blockStart, i, byteOffset, bitOffset,
						// cramRecord.getAlignmentStart());

					} else
						cramRecord = reader.readRecord();

				} catch (Exception e) {
					log.error("Failed to read record: " + i);
					if (cramRecord != null)
						log.error(cramRecord.toString());
					log.error(e);
					throw e;
				}
			}

			// if (blockCounter++ > 2)
			// break;
			long blockTime2 = System.currentTimeMillis();
			log.debug(String.format("Block read in %.2f seconds.", (blockTime2 - blockTime1) / 1000f));
		}

		indexPS.flush();
		long time2 = System.currentTimeMillis();
		log.info("Decoded in: " + (time2 - time1) + " millis");
	}

	@Parameters(commandDescription = "Build CRAM index. ")
	static class Params {
		@Parameter(names = { "--input-cram-file" }, converter = FileConverter.class, description = "The path to the CRAM file to be indexed.")
		File cramFile;

		@Parameter(names = { "--max-sequences" }, description = "Stop after indexing this many reference sequences (chromosomes).")
		int maxSequences = Integer.MAX_VALUE;

		@Parameter(names = { "--max-records" }, description = "Stop after indexing this many records.")
		long maxRecords = Long.MAX_VALUE;

		@Parameter(names = { "--reference-fasta-file" }, converter = FileConverter.class, description = "The path to the reference fasta file, uncompressed and indexed.")
		File reference;

		@Parameter(names = { "--index-file" }, converter = FileConverter.class, description = "The path for the output index file. By default a file with '.crai' extention will be created.")
		File indexFile;

		@Parameter(names = { "-h", "--help" }, description = "Print help and quit")
		boolean help = false;

		@Parameter(names = { "--resolution" }, description = "Index one record per this many. Bigger value means lower resolution.")
		int resolution = 1000;

	}
}
