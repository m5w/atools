h_sources = alignment_io.h \
	    array2d.h \
	    exception.h \
	    exception_type.h \
	    filelib.h \
	    getopt_long_exception.h \
	    gzstream.h \
	    null_deleter.h \
	    optional.h
cc_sources = alignment_io.cc \
	     exception_type.cc \
	     filelib.cc \
	     gzstream.cc
bin_PROGRAMS = atools
lib_LIBRARIES = libutils.a
libutils_a_SOURCES = $(h_sources) $(cc_sources)
atools_SOURCES = atools.cc
atools_LDADD = libutils.a $(ZLIB_LIBS)
