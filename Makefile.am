NULL =
ACLOCAL_AMFLAGS = -I m4
bin_PROGRAMS = pandaseq-sam
lib_LTLIBRARIES = libpandaseq-sam.la
library_includedir=$(includedir)/$(LIB_NAME)
library_include_HEADERS = pandaseq-sam.h
if PTHREAD
library_include_HEADERS += pandaseq-sam-mux.h
endif
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = $(LIB_NAME).pc
vapidir = $(datadir)/vala/vapi
dist_vapi_DATA = $(LIB_NAME).vapi $(LIB_NAME).deps
man1_MANS = pandaseq-sam.1
docdir = $(datadir)/doc/@PACKAGE@
doc_DATA = README
pandaseq_sam_CPPFLAGS = $(PANDASEQ_CFLAGS) -std=c99 -pedantic -g
pandaseq_sam_LDFLAGS = $(PANDASEQ_LIBS)
pandaseq_sam_LDADD = libpandaseq-sam.la
pandaseq_sam_SOURCES = main.c
libpandaseq_sam_la_CPPFLAGS = \
	$(HTS_CFLAGS) \
	$(PANDASEQ_CFLAGS) \
	$(PTHREAD_CFLAGS) \
	-DPANDA_SAM_LIB_COMPILING \
	-g \
	-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE \
	$(NULL)
libpandaseq_sam_la_LDFLAGS = \
	$(HTS_LIBS) \
	$(PANDASEQ_LIBS) \
	$(PTHREAD_LIBS) \
	-export-symbols-regex '^panda_' \
 	-version-info $(LIB_VER) \
	-no-undefined \
	$(NULL)
libpandaseq_sam_la_SOURCES = \
	args.c \
	reader.c \
	seqid.c \
	support.c \
	$(NULL)

