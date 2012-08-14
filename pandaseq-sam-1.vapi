/**
 * PANDAseq Assembler for Paired SAM/BAM Illumina
 */
[CCode(cheader_filename = "pandaseq-sam.h")]
namespace Panda.SAM {
	/**
	 * Get the version of the underlying SAM tools library used.
	 */
	[CCode(cname = "panda_sam_version")]
	public unowned string version();

	/**
	 * Parse a SAM/BAM style Illumina sequence.
	 *
	 * These are different from normal Illumina names as the tags have been stripped.
	 */
	[CCode(cname = "panda_seqid_parse_sam")]
	public bool parse_id(out identifier id, string str);

	/**
	 * Create an object to read sequences from a SAM file
	 *
	 * @param filename the filename containing paired-end Illumina sequences
	 * @param logger the logging function to use during assembly.
	 * @param binary whether the file is binary (BAM) or text (SAM)
	 * @param tag a tag to replace the missing Illumina barcoding tag
	 */
	[CCode(cname = "panda_create_sam_reader")]
	public NextSeq? create_reader(string filename, Logger logger, bool binary, string? tag = null);
	/**
	 * Create a new assembler for given a SAM file.
	 * @see create_reader
	 */
	[CCode(cname = "panda_assembler_open_sam")]
	public Assembler? open(string filename, owned Logger logger, bool binary, string? tag = null);
	/**
	 * Create a new multiplexed reader for given a SAM file.
	 * @see create_reader
	 */
	[CCode(cname = "panda_mux_open_sam")]
	public Mux? open_mux(string filename, owned Logger logger, bool binary, string? tag = null);
}
