AUTOMAKE_OPTIONS = foreign

#  Header files
includedir= @includedir@/cteq

include_HEADERS = cteqpdf.h

#  Core library : libcteqpdf.la
lib_LTLIBRARIES = libcteqpdf.la 

libcteqpdf_la_SOURCES = cteqpdf.c cteq-database.c cteqpdf-f77.c
libcteqpdf_la_CPPFLAGS = -DCTEQ_TBL_PATH=\"$(datadir)/cteq-pdf\"
libcteqpdf_la_LDFLAGS = -release $(VERSION)
libcteqpdf_la_LIBADD = -lm 

