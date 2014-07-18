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

#include "config.h"
#include <errno.h>
#include <unistd.h>

#include "pandaseq-sam.h"
#include <htslib/hts.h>
#include <htslib/khash.h>
#include <htslib/sam.h>

KHASH_MAP_INIT_STR(seq, bam1_t *)

struct reader_data {
	htsFile *file;
	 khash_t(
		seq) * pool;
	PandaLogProxy logger;
	panda_qual *forward;
	size_t forward_length;
	panda_qual *reverse;
	size_t reverse_length;
	size_t tag_length;
	char tag[PANDA_TAG_LEN];
	bam_hdr_t *header;
	FILE *orphan_file;
};

void ps_fill(
	bam1_t *bam,
	panda_qual *seq,
	size_t *seq_length) {
	size_t it;
	*seq_length = bam->core.l_qseq;
	for (it = 0; it < *seq_length; it++) {
		size_t pos = (bam->core.flag & BAM_FREVERSE) ? (*seq_length - it - 1) : it;
		seq[pos].nt = (panda_nt) (bam_seqi(bam_get_seq(bam), it));
		seq[pos].qual = bam_get_qual(bam)[it];
	}
}

bool damaged_seq(
	bam1_t *seq,
	PandaCode *code) {
	if (seq->core.l_qseq < 1) {
		*code = PANDA_CODE_NO_DATA;
		return true;
	}
	if ((size_t) seq->core.l_qseq > PANDA_MAX_LEN) {
		*code = PANDA_CODE_READ_TOO_LONG;
		return true;
	}
	if (!(seq->core.flag & BAM_FPAIRED)) {
		*code = PANDA_CODE_NOT_PAIRED;
		return true;
	}
	return false;
}

#define show_flag(x) if (seq->core.flag & x) fprintf(data->orphan_file, " " #x);

void write_orphan(
	struct reader_data *data,
	bam1_t *seq,
	PandaCode seq_err) {
	if (data->orphan_file != NULL) {
		size_t it;
		fprintf(data->orphan_file, "@%s", bam_get_qname(seq));
		show_flag(BAM_FPAIRED);
		show_flag(BAM_FPROPER_PAIR);
		show_flag(BAM_FUNMAP);
		show_flag(BAM_FMUNMAP);
		show_flag(BAM_FREVERSE);
		show_flag(BAM_FMREVERSE);
		show_flag(BAM_FREAD1);
		show_flag(BAM_FREAD2);
		show_flag(BAM_FSECONDARY);
		show_flag(BAM_FQCFAIL);
		show_flag(BAM_FDUP);
		show_flag(BAM_FSUPPLEMENTARY);
		fprintf(data->orphan_file, "\n");
		for (it = 0; it < (size_t) seq->core.l_qseq; it++) {
			fputc(panda_nt_to_ascii((panda_nt) (bam_seqi(bam_get_seq(seq), it))), data->orphan_file);
		}
		fprintf(data->orphan_file, "\n+\n");
		for (it = 0; it < (size_t) seq->core.l_qseq; it++) {
			fputc(33 + bam_get_qual(seq)[it], data->orphan_file);
		}
		fprintf(data->orphan_file, "\n");
	} else if (panda_debug_flags & PANDA_DEBUG_FILE) {
		panda_log_proxy_write(data->logger, seq_err, NULL, NULL, bam_get_qname(seq));
	}
}

bool ps_next(
	panda_seq_identifier *id,
	panda_qual **forward,
	size_t *forward_length,
	panda_qual **reverse,
	size_t *reverse_length,
	struct reader_data *data) {
	int res;
	khiter_t key;
	bam1_t *seq = bam_init1();

	*forward = NULL;
	*forward_length = 0;
	*reverse = NULL;
	*reverse_length = 0;
	while ((res = sam_read1(data->file, data->header, seq)) >= 0) {
		PandaCode seq_err;
		if (damaged_seq(seq, &seq_err)) {
			write_orphan(data, seq, seq_err);
			continue;
		}
		key = kh_get(seq, data->pool, bam_get_qname(seq));
		if (key == kh_end(data->pool)) {
			int ret;
			khiter_t key;
			key = kh_put(seq, data->pool, bam_get_qname(seq), &ret);
			if (ret == 0) {
				if (panda_debug_flags & PANDA_DEBUG_FILE) {
					panda_log_proxy_write(data->logger, PANDA_CODE_PREMATURE_EOF, NULL, NULL, bam_get_qname(seq));
				}
				bam_destroy1(seq);
				return false;
			}
			kh_value(data->pool, key) = seq;
			seq = bam_init1();
		} else {
			bam1_t *mate = kh_value(data->pool, key);
			kh_del(seq, data->pool, key);

			if (!panda_seqid_parse_sam(id, bam_get_qname(seq))) {
				if (panda_debug_flags & PANDA_DEBUG_FILE) {
					panda_log_proxy_write(data->logger, PANDA_CODE_ID_PARSE_FAILURE, NULL, NULL, bam_get_qname(seq));
				}
				bam_destroy1(mate);
				bam_destroy1(seq);
				return false;
			}
			memcpy(id->tag, data->tag, data->tag_length + 1);

			if (seq->core.flag & BAM_FREAD1) {
				ps_fill(seq, data->forward, &data->forward_length);
				ps_fill(mate, data->reverse, &data->reverse_length);
			} else {
				ps_fill(mate, data->forward, &data->forward_length);
				ps_fill(seq, data->reverse, &data->reverse_length);
			}

			bam_destroy1(seq);
			bam_destroy1(mate);
			*forward = data->forward;
			*forward_length = data->forward_length;
			*reverse = data->reverse;
			*reverse_length = data->reverse_length;
			return true;
		}
	}
	/* -1 is normal end of file. */
	if (res < -1 && panda_debug_flags & PANDA_DEBUG_FILE) {
		panda_log_proxy_write(data->logger, PANDA_CODE_PREMATURE_EOF, NULL, NULL, bam_get_qname(seq));
	}
	bam_destroy1(seq);
	return false;
}

void ps_destroy(
	struct reader_data *data) {
	khiter_t key;
	bam_hdr_destroy(data->header);
	hts_close(data->file);
	for (key = kh_begin(data->pool); key != kh_end(data->pool); key++) {
		if (kh_exist(data->pool, key)) {
			bam1_t *seq = kh_value(data->pool, key);
			write_orphan(data, seq, PANDA_CODE_PARSE_FAILURE);
			bam_destroy1(seq);
		}
	}
	kh_destroy(seq, data->pool);
	panda_log_proxy_unref(data->logger);
	if (data->orphan_file != NULL) {
		fclose(data->orphan_file);
	}
	free(data->forward);
	free(data->reverse);
	free(data);
}

PandaNextSeq panda_create_sam_reader(
	const char *filename,
	PandaLogProxy logger,
	bool binary,
	const char *tag,
	void **user_data,
	PandaDestroy *destroy) {
	return panda_create_sam_reader_orphans(filename, logger, binary, tag, NULL, user_data, destroy);
}

PandaNextSeq panda_create_sam_reader_orphans(
	const char *filename,
	PandaLogProxy logger,
	bool binary,
	const char *tag,
	const char *orphan_file,
	void **user_data,
	PandaDestroy *destroy) {
	struct reader_data *data;

	*destroy = NULL;
	*user_data = NULL;

	data = malloc(sizeof(struct reader_data));
	if (data == NULL) {
		return NULL;
	}
	data->forward = calloc(PANDA_MAX_LEN, sizeof(panda_qual));
	if (data->forward == NULL) {
		free(data);
		return NULL;
	}

	data->reverse = calloc(PANDA_MAX_LEN, sizeof(panda_qual));
	if (data->reverse == NULL) {
		free(data->forward);
		free(data);
		return NULL;
	}

	data->file = hts_open(filename, binary ? "rb" : "r");
	if (data->file == NULL) {
		free(data->forward);
		free(data->reverse);
		free(data);
		return NULL;
	}
	if (tag == NULL) {
		data->tag_length = 0;
		data->tag[0] = '\0';
	} else {
		data->tag_length = strlen(tag);
		if (data->tag_length >= PANDA_TAG_LEN) {
			data->tag_length = PANDA_TAG_LEN - 1;
		}
		memcpy(data->tag, tag, data->tag_length);
		data->tag[data->tag_length] = '\0';
	}
	if (orphan_file == NULL) {
		data->orphan_file = NULL;
	} else {
		data->orphan_file = NULL;
		if (access(orphan_file, F_OK) == -1 && errno == ENOENT) {
			data->orphan_file = fopen(orphan_file, "w");
			if (data->orphan_file == NULL) {
				perror(orphan_file);
			}
		} else {
			fprintf(stderr, "%s: refusing to overwrite file.\n", orphan_file);
		}
		if (data->orphan_file == NULL) {
			hts_close(data->file);
			free(data->forward);
			free(data->reverse);
			free(data);
			return NULL;
		}
	}
	data->pool = kh_init(seq);
	data->header = sam_hdr_read(data->file);
	data->logger = panda_log_proxy_ref(logger);
	*destroy = (PandaDestroy) ps_destroy;
	*user_data = data;
	return (PandaNextSeq) ps_next;
}
