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
#define _PANDASEQ_SAM_H
#include <pandaseq.h>

/**
 * Get the version of the underlying SAM tools library used.
 */
const char *panda_sam_version(void);

/**
 * Parse a SAM/BAM style Illumina sequence.
 *
 * These are different from normal Illumina names as the tags have been stripped.
 */
bool panda_seqid_parse_sam(panda_seq_identifier *id, char *str);

/**
 * Create an object to read sequences from a SAM file
 *
 * @param filename the filename containing paired-end Illumina sequences
 * @param logger, logger_data the logging function to use during assembly. The logging function will not be memory managed.
 * @param binary whether the file is binary (BAM) or text (SAM)
 * @param qualmin the value to strip from the quality scores. Usually 33 or 64, depending on CASAVA version.
 * @param user_data where to store the user_data for this function
 * @param destroy where to store the destroy function for the user data
 */
PandaNextSeq panda_create_sam_reader(/*@notnull@ */ char *filename,
				       /*@notnull@ */ PandaLogger logger,
				       /*@null@ */ void *logger_data,
				       bool binary,
				       unsigned char qualmin,
				       /*@notnull@@out@ */ void **user_data,
				       /*@notnull@@out@ */
				       PandaDestroy * destroy);
/**
 * Create a new assembler for given a SAM file.
 * @see panda_create_sam_reader
 */
PandaAssembler panda_assembler_open_sam( /*@notnull@ */ char *filename,
								     /*@notnull@ */ PandaLogger logger,
								     /*@null@ */ void *logger_data,
								     /*@null@ */ PandaDestroy logger_destroy,
								     bool binary,
								     unsigned char qualmin);
/**
 * Create a new multiplexed reader for given a SAM file.
 * @see panda_create_sam_reader
 */
PandaMux panda_mux_open_sam( /*@notnull@ */ char *filename,
						  /*@notnull@ */ PandaLogger logger,
						  /*@null@ */ void *logger_data,
						  /*@null@ */ PandaDestroy logger_destroy,
						  bool binary,
						  unsigned char qualmin);
#endif
