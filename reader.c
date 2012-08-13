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

#include "pandaseq-sam.h"
#include <khash.h>
#include <sam.h>
KHASH_MAP_INIT_STR(seq, bam1_t*)

struct reader_data {
	samfile_t *file;
	khash_t(seq) *pool;
	PandaLogger logger;
	void *logger_data;
	unsigned char qualmin;
	size_t forward_length;
	size_t reverse_length;
	panda_qual forward[PANDA_MAX_LEN];
	panda_qual reverse[PANDA_MAX_LEN];
	char tag[PANDA_TAG_LEN];
	size_t tag_length;
};

bool ps_next(panda_seq_identifier *id, panda_qual **forward, size_t *forward_length, panda_qual **reverse, size_t *reverse_length, struct reader_data *data) {
	size_t it;
	int res;
	khiter_t key;
	bam1_t *seq = bam_init1();
	while ((res = samread(data->file, seq)) > 0) {
		if (seq->core.l_qseq < 1 || seq->core.l_qseq > PANDA_MAX_LEN || !(seq->core.flag & BAM_FPAIRED)) {
			continue;
		}
		key = kh_get(seq, data->pool, bam1_qname(seq));
		if (key == kh_end(data->pool)) {
			int ret;
			khiter_t key;
			key = kh_put(seq, data->pool, bam1_qname(seq), &ret);
			if (ret == 0) {
				if (data->logger != NULL && panda_debug_flags & PANDA_DEBUG_FILE) {
					data->logger(PANDA_CODE_PREMATURE_EOF, NULL, bam1_qname(seq), data->logger_data);
				}
				bam_destroy1(seq);
				return false;
			}
			kh_value(data->pool, key) = seq;
			seq = bam_init1();
		} else {
			bam1_t *bam_forward;
			bam1_t *bam_reverse;
			bam1_t *mate = kh_value(data->pool, key);
			kh_del(seq, data->pool, key);

			if (!panda_seqid_parse_sam(id, bam1_qname(seq))) {
				bam_destroy1(mate);
				bam_destroy1(seq);
				if (data->logger != NULL && panda_debug_flags & PANDA_DEBUG_FILE) {
					data->logger(PANDA_CODE_ID_PARSE_FAILURE, NULL, bam1_qname(seq), data->logger_data);
				}
				return false;
			}
			memcpy(id->tag, data->tag, data->tag_length + 1);

			if (seq->core.flag & BAM_FREAD1) {
				bam_forward = seq;
				bam_reverse = mate;
			} else {
				bam_forward = mate;
				bam_reverse = seq;
			}

			data->forward_length = bam_forward->core.l_qseq;
			for(it = 0; it < data->forward_length; it++) {
				data->forward[it].nt = (panda_nt)(bam1_seqi(bam1_seq(bam_forward), it));
				data->forward[it].qual = bam1_qual(bam_forward)[it] - data->qualmin;
			}

			data->reverse_length = bam_reverse->core.l_qseq;
			for(it = 0; it < data->reverse_length; it++) {
				data->reverse[it].nt = (panda_nt)bam1_seqi(bam1_seq(bam_reverse), it);
				data->reverse[it].qual = bam1_qual(bam_reverse)[it] - data->qualmin;
			}
			bam_destroy1(seq);
			bam_destroy1(mate);
			return true;
		}
	}
	if (res == -2 && data->logger != NULL && panda_debug_flags & PANDA_DEBUG_FILE) {
		data->logger(PANDA_CODE_PREMATURE_EOF, NULL, bam1_qname(seq), data->logger_data);
	}
	bam_destroy1(seq);
	return false;
}

void ps_destroy(struct reader_data *data) {
	khiter_t key;
	samclose(data->file);
	for (key = kh_begin(data->pool); key != kh_end(data->pool); key++) {
		if (kh_exist(data->pool, key)) {
			bam_destroy1(kh_value(data->pool, key));
		}
	}
	kh_destroy(seq, data->pool);
	free(data);
}

PandaNextSeq panda_create_sam_reader(char *filename, PandaLogger logger, void *logger_data, bool binary, char *tag, unsigned char qualmin, void **user_data, PandaDestroy *destroy) {
	struct reader_data *data;

	*destroy = NULL;
	*user_data = NULL;

	data = malloc(sizeof(struct reader_data));
	if (data == NULL) {
		return NULL;
	}

	data->file = samopen(filename, binary ? "rb" : "r", NULL);
	if (data->file == NULL) {
		free(data);
		return NULL;
	}
	data->qualmin = qualmin;
	data->pool = kh_init(seq);
	data->qualmin = qualmin;
	if (tag == NULL) {
		data->tag[0] = '\0';
	} else {
		data->tag_length = strlen(tag);
		if (data->tag_length >= PANDA_TAG_LEN) {
			data->tag_length = PANDA_TAG_LEN - 1;
		}
		memcpy(data->tag, tag, data->tag_length);
		data->tag[data->tag_length] = '\0';
	}
	*destroy = (PandaDestroy)ps_destroy;
	*user_data = data;
	return (PandaNextSeq) ps_next;
}
