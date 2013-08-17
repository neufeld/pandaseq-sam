/* PANDAseq -- Assemble paired SAM/BAM Illumina reads and strip the region between amplification primers.
     Copyright (C) 2012  Andre Masella

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef _PANDASEQ_SAM_H
#        define _PANDASEQ_SAM_H
#        include <pandaseq.h>
#        ifdef __cplusplus
#                define EXTERN_C_BEGIN  extern "C" {
#                define EXTERN_C_END    }
#        else
#                define EXTERN_C_BEGIN
#                define EXTERN_C_END
#        endif
#        ifdef _WIN32
#                ifdef PANDA_SAM_LIB_COMPILING
#                        define PANDA_SAM_EXTERN extern __declspec(dllexport)
#                else
#                        define PANDA_SAM_EXTERN extern __declspec(dllimport)
#                endif
#        else
#                define PANDA_SAM_EXTERN extern
#        endif
EXTERN_C_BEGIN
/**
 * Get the version of the underlying SAM tools library used.
 */
const char *panda_sam_version(
	void);

/**
 * Parse a SAM/BAM style Illumina sequence.
 *
 * These are different from normal Illumina names as the tags have been stripped.
 */
bool panda_seqid_parse_sam(
	panda_seq_identifier *id,
	char *str);

/**
 * Create an object to read sequences from a SAM file
 *
 * @filename: the filename containing paired-end Illumina sequences
 * @logger: the logging to use during assembly
 * @binary: whether the file is binary (BAM) or text (SAM)
 * @tag:(allow-none): a tag to replace the missing Illumina barcoding tag
 * Returns:(closure user_data) (scope notified): a sequence source callback
 */
PandaNextSeq panda_create_sam_reader(
	const char *filename,
	PandaLogProxy logger,
	bool binary,
	const char *tag,
	void **user_data,
	PandaDestroy *destroy);
/**
 * Create a new assembler for given a SAM file.
 * @see panda_create_sam_reader
 */
PandaAssembler panda_assembler_open_sam(
	const char *filename,
	PandaLogProxy logger,
	bool binary,
	const char *tag);
/**
 * Create a new multiplexed reader for given a SAM file.
 * @see panda_create_sam_reader
 */
PandaMux panda_mux_open_sam(
	const char *filename,
	PandaLogProxy logger,
	bool binary,
	const char *tag);
/**
 * The standard argument handler for a SAM file of pair-end reads.
 */
typedef struct panda_args_sam *PandaArgsSam;

/**
 * Command line arguments for a SAM file of pair-end reads.
 */
PANDA_SAM_EXTERN const panda_tweak_general *const panda_args_sam_args[];
PANDA_SAM_EXTERN const size_t panda_args_sam_args_length;

/**
 * Create a new argument handler.
 */
PandaArgsSam panda_args_sam_new(
	);

/**
 * Cleanup the argument handler.
 */
void panda_args_sam_free(
	PandaArgsSam data);

/**
 * Process the command line arguments for the SAM argument handler.
 */
bool panda_args_sam_tweak(
	PandaArgsSam data,
	char flag,
	const char *argument);

/**
 * Initialise the sequence stream for the SAM argument handler.
 */
PandaNextSeq panda_args_sam_opener(
	PandaArgsSam data,
	PandaLogProxy logger,
	PandaFailAlign *fail,
	void **fail_data,
	PandaDestroy *fail_destroy,
	void **next_data,
	PandaDestroy *next_destroy);

/**
 * Do additional assembly setup for the SAM argument handler.
 */
bool panda_args_sam_setup(
	PandaArgsSam data,
	PandaAssembler assembler);

EXTERN_C_END
#endif
