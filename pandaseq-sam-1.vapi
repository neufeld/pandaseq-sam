/**
 * PANDAseq Assembler for Paired SAM/BAM Illumina
 */
[CCode (cheader_filename = "pandaseq-sam.h")]
namespace Panda.SAM {
	/**
	 * Get the version of the underlying SAM tools library used.
	 */
	[CCode (cname = "panda_sam_version")]
	public unowned string version ();

	/**
	 * Parse a SAM/BAM style Illumina sequence.
	 *
	 * These are different from normal Illumina names as the tags have been stripped.
	 */
	[CCode (cname = "panda_seqid_parse_sam")]
	public bool parse_id (out identifier id, string str);

	/**
	 * Create an object to read sequences from a SAM file
	 *
	 * @param filename the filename containing paired-end Illumina sequences
	 * @param logger the logging to use during assembly.
	 * @param binary whether the file is binary (BAM) or text (SAM)
	 * @param tag a tag to replace the missing Illumina barcoding tag
	 * @param orphan_file the FASTQ file where unpaired/damaged/broken reads should be place.
	 * @param reverse_direction reverse the direction of the reverse read. Normally true.
	 */
	[CCode (cname = "panda_create_sam_reader_ex")]
	public NextSeq? create_reader (string filename, LogProxy logger, bool binary, string? tag = null, string? orphans_file = null, bool reverse_direction = true);
	/**
	 * Create a new assembler for given a SAM file.
	 * @see create_reader
	 */
	[CCode (cname = "panda_assembler_open_sam")]
	public Assembler? open (string filename, LogProxy logger, bool binary, string? tag = null);
	/**
	 * Create a new multiplexed reader for given a SAM file.
	 * @see create_reader
	 */
	[CCode (cname = "panda_mux_open_sam", cheader_filename = "pandaseq-sam-mux.h")]
	public Mux? open_mux (string filename, LogProxy logger, bool binary, string? tag = null);

	/**
	 * The standard argument handler for a SAM file of pair-end reads.
	 */
	[CCode (cname = "struct panda_args_sam", free_function = "panda_args_sam_free")]
	[Compact]
	public class Args {

		/**
		 * Command line arguments for a SAM file of pair-end reads.
		 */
		[CCode (cname = "panda_args_sam_args", array_length_cexpr = "panda_args_sam_args_length", array_length_type = "size_t")]
		extern const Tweak.general? [] args;

		/**
		 * Create a new argument handler.
		 */
		[CCode (cname = "panda_args_sam_new")]
		public Args ();

		/**
		 * Process the command line arguments for the SAM argument handler.
		 */
		[CCode (cname = "panda_args_sam_tweak")]
		public bool tweak (char flag, string argument);

		/**
		 * Initialise the sequence stream for the SAM argument handler.
		 */
		[CCode (cname = "panda_args_sam_opener")]
		public NextSeq opener (LogProxy logger, out FailAlign? fail);

		/**
		 * Do additional assembly setup for the SAM argument handler.
		 */
		[CCode (cname = "panda_args_sam_setup")]
		public bool setup (Assembler assembler);
	}
}
