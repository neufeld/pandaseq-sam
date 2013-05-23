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
#include<stdio.h>
#include<stdlib.h>
#include "config.h"
#ifdef HAVE_PTHREAD
#        include<pthread.h>
#endif
#include "pandaseq-sam.h"

int main(
	int argc,
	char **argv) {
	PandaOutputSeq output;
	void *output_data;
	PandaDestroy output_destroy;
	PandaAssembler assembler;
	PandaArgsSam data = panda_args_sam_new();
	PandaMux mux;
	bool result;
	int threads;

	if (!panda_parse_args(argv, argc, panda_stdargs, panda_stdargs_length, panda_args_sam_args, panda_args_sam_args_length, (PandaTweakGeneral) panda_args_sam_tweak, (PandaOpener) panda_args_sam_opener, (PandaSetup) panda_args_sam_setup, data, &assembler, &mux, &threads, &output, &output_data, &output_destroy)) {
		panda_args_sam_free(data);
		return 1;
	}
	result = panda_run_pool(threads, assembler, mux, output, output_data, output_destroy);
	panda_args_sam_free(data);
	return result ? 0 : 1;
}
