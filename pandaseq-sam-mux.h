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

#ifndef _PANDASEQ_SAM_MUX_H
#        define _PANDASEQ_SAM_MUX_H
#        include <pandaseq-sam.h>
#        include <pandaseq-mux.h>
#        ifdef __cplusplus
#                define EXTERN_C_BEGIN  extern "C" {
#                define EXTERN_C_END    }
#        else
#                define EXTERN_C_BEGIN
#                define EXTERN_C_END
#        endif
EXTERN_C_BEGIN
/**
 * Create a new multiplexed reader for given a SAM file.
 * @see panda_create_sam_reader
 */
PandaMux panda_mux_open_sam(
	const char *filename,
	PandaLogProxy logger,
	bool binary,
	const char *tag);
EXTERN_C_END
#endif
