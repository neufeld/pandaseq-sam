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
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif
#include <sam.h>

const char *panda_sam_version(void) {
	return BAM_VERSION;
}

PandaAssembler panda_assembler_open_sam(char *filename, PandaLogger logger, void *logger_data, PandaDestroy logger_destroy, bool binary, unsigned char qualmin){
	void *next_data;
	PandaDestroy next_destroy;
	PandaNextSeq next;
	next = panda_create_sam_reader(filename, logger, logger_data, binary, qualmin, &next_data, &next_destroy);
	if (next == NULL) {
		return NULL;
	}
	return panda_assembler_new(next, next_data, next_destroy, logger, logger_data, logger_destroy);
}

#ifdef HAVE_PTHREAD
PandaMux panda_mux_open_sam(char *filename, PandaLogger logger, void *logger_data, PandaDestroy logger_destroy, bool binary, unsigned char qualmin) {
	void *next_data;
	PandaDestroy next_destroy;
	PandaNextSeq next;
	next = panda_create_sam_reader(filename, logger, logger_data, binary, qualmin, &next_data, &next_destroy);
	if (next == NULL) {
		return NULL;
	}
	return panda_mux_new(next, next_data, next_destroy, logger, logger_data, logger_destroy);
}
#endif
