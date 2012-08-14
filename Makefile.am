NULL =
ACLOCAL_AMFLAGS = -I m4
bin_PROGRAMS = pandaseq-sam
lib_LTLIBRARIES = libpandaseq-sam.la
library_includedir=$(includedir)/$(LIB_NAME)
library_include_HEADERS = pandaseq-sam.h
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = $(LIB_NAME).pc
vapidir = $(datadir)/vala/vapi
dist_vapi_DATA = $(LIB_NAME).vapi $(LIB_NAME).deps
man1_MANS = pandaseq-sam.1
docdir = $(datadir)/doc/@PACKAGE@
doc_DATA = README
pandaseq_sam_CPPFLAGS = $(PANDASEQ_CFLAGS) -std=gnu99 -pedantic -g
pandaseq_sam_LDFLAGS = $(PANDASEQ_LIBS) libpandaseq-sam.la
pandaseq_sam_SOURCES = main.c
libpandaseq_sam_la_CPPFLAGS = \
	$(LIB_LIBS) \
	$(PANDASEQ_CFLAGS) \
	$(PTHREAD_CFLAGS) \
	-g \
	-Isamtools \
	-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE \
	$(NULL)
libpandaseq_sam_la_CFLAGS = $(LIB_CFLAGS)
libpandaseq_sam_la_LDFLAGS = \
	$(LIB_CPPFLAGS) \
	$(PANDASEQ_LIBS) \
	$(PTHREAD_LIBS) \
	-export-symbols-regex '^panda_' \
	$(NULL)
libpandaseq_sam_la_SOURCES = \
	reader.c \
	seqid.c \
	support.c \
	samtools/bam_aux.c \
	samtools/bam.c \
	samtools/bam_cat.c \
	samtools/bam_import.c \
	samtools/bam_index.c \
	samtools/bam_lpileup.c \
	samtools/bam_md.c \
	samtools/bam_pileup.c \
	samtools/bam_reheader.c \
	samtools/bam_sort.c \
	samtools/bedidx.c \
	samtools/bgzf.c \
	samtools/faidx.c \
	samtools/knetfile.c \
	samtools/kprobaln.c \
	samtools/kstring.c \
	samtools/razf.c \
	samtools/sam.c \
	samtools/sam_header.c \
	$(NULL)

