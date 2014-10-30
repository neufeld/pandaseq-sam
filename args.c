/* PANDAseq-SAM -- Assemble paired SAM-format Illumina reads and strip the region between amplification primers.
     Copyright (C) 2011-2013  Andre Masella

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
#define _POSIX_C_SOURCE 2
#include<ctype.h>
#include<errno.h>
#include<float.h>
#include<limits.h>
#include<stdbool.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<unistd.h>
#include "config.h"
#include "pandaseq-sam.h"
#ifdef HAVE_PTHREAD
#        include<pthread.h>
#endif

struct panda_args_sam {
	bool binary;
	const char *filename;
	const char *orphans_file;
	char tag[PANDA_TAG_LEN];
	bool no_algn_qual;
	PandaWriter no_algn_writer;
};

PandaArgsSam panda_args_sam_new(
	) {
	PandaArgsSam data = malloc(sizeof(struct panda_args_sam));
	data->binary = false;
	data->filename = NULL;
	data->orphans_file = NULL;
	data->tag[0] = '\0';
	data->no_algn_qual = false;
	data->no_algn_writer = NULL;
	return data;
}

void panda_args_sam_free(
	PandaArgsSam data) {
	panda_writer_unref(data->no_algn_writer);
	free(data);
}

bool panda_args_sam_tweak(
	PandaArgsSam data,
	char flag,
	const char *argument) {
	switch (flag) {
	case 'b':
		data->binary = true;
		return true;
	case 'B':
		strncpy(data->tag, argument, PANDA_TAG_LEN);
		if (data->tag[PANDA_TAG_LEN - 1] != '\0') {
			fprintf(stderr, "Replacement tag %s is too long.", argument);
			return false;
		}
		return true;
	case 'r':
		data->orphans_file = argument;
		return true;
	case 'f':
		data->filename = (strcmp(optarg, "-") == 0) ? "/dev/stdin" : argument;
		return true;
	case 'u':
	case 'U':
		data->no_algn_qual = flag == 'U';
		panda_writer_unref(data->no_algn_writer);
		data->no_algn_writer = panda_writer_open_file(argument, false);
		if (data->no_algn_writer == NULL)
			perror(argument);
		return (data->no_algn_writer != NULL);
	default:
		return false;
	}
}

#define MAYBE(x) if (x != NULL) *x

PandaNextSeq panda_args_sam_opener(
	PandaArgsSam data,
	PandaLogProxy logger,
	PandaFailAlign *fail,
	void **fail_data,
	PandaDestroy *fail_destroy,
	void **next_data,
	PandaDestroy *next_destroy) {

	if (data->no_algn_writer != NULL) {
		*fail = (PandaFailAlign) (data->no_algn_qual ? panda_output_fail_qual : panda_output_fail);
		*fail_data = data->no_algn_writer;
		*fail_destroy = (PandaDestroy) panda_writer_unref;
		data->no_algn_writer = NULL;
	} else {
		*fail = NULL;
		*fail_data = NULL;
		*fail_destroy = NULL;
	}

	if (data->filename == NULL) {
		MAYBE(next_data) = NULL;
		MAYBE(next_destroy) = NULL;
		return false;
	}
	return panda_create_sam_reader_ex(data->filename, logger, data->binary, data->tag, data->orphans_file, next_data, next_destroy);
}

bool panda_args_sam_setup(
	PandaArgsSam data,
	PandaAssembler assembler) {

	(void) data;
	(void) assembler;
	return true;
}

const panda_tweak_general args_filename = { 'f', false, "file.sam", "Input SAM/BAM file containing forward reads.", false };

const panda_tweak_general args_bin = { 'b', true, NULL, "Read a binary (BAM) file rather than a text (SAM) file.", false };

static const panda_tweak_general args_unalign_qual = { 'U', true, "unaligned.txt", "File to write unalignable read pairs with quality scores.", false };

const panda_tweak_general args_code = { 'B', true, "code", "Replace the Illumina multiplexing barcode stripped during processing into SAM/BAM.", false };

const panda_tweak_general args_orphans = { 'r', true, "orphans.fastq", "Write all reads from the SAM/BAM that could not be paired or were discarded to a FASTQ file.", false };

static const panda_tweak_general args_unalign = { 'u', true, "unaligned.txt", "File to write unalignable read pairs.", false };

const panda_tweak_general *const panda_args_sam_args[] = {
	&args_code,
	&args_bin,
	&args_unalign_qual,
	&args_filename,
	&args_orphans,
	&args_unalign
};

const size_t panda_args_sam_args_length = sizeof(panda_args_sam_args) / sizeof(panda_tweak_general *);
