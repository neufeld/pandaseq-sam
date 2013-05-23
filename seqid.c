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
#include <string.h>

#include "pandaseq-sam.h"

#define PARSE_CHUNK if (*input == '\0') return 0; for(;*input != '\0' && *input != ':' && *input != ' '; input++)
#define PARSE_INT do { value = 0; PARSE_CHUNK { if (*input >= '0' && *input <= '9') { value = 10*value + (*input - '0'); } else { return 0; } } } while(0)

bool panda_seqid_parse_sam(
	panda_seq_identifier *id,
	char *input) {
	char *dest;
	int value;
	if (strpbrk(input, "/#") != NULL) {
		return false;
	}
	id->run = 0;
	id->flowcell[0] = '\0';
	id->tag[0] = '\0';
	dest = id->instrument;
	PARSE_CHUNK {
		*dest++ = (*input);
	}
	input++;
	*dest = '\0';
	PARSE_INT;
	input++;
	id->lane = value;
	PARSE_INT;
	input++;
	id->tile = value;
	PARSE_INT;
	input++;
	id->x = value;
	PARSE_INT;
	id->y = value;
	return *input == '\0';
}
