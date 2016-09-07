#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.63)
AC_INIT([cteq-pdf],[1.0.4],[Zoltan.Nagy@cern.ch])
AM_INIT_AUTOMAKE(cteq-pdf, 1.0.4)

AC_CONFIG_SRCDIR([src/cteqpdf.c])
AC_CONFIG_HEADER([config.h])
AC_PREFIX_DEFAULT(${HOME})

# Checks for programs.
LT_INIT
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m], [cos])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow])


AC_CONFIG_FILES([Makefile src/Makefile
cteq-tbls/Makefile
cteq-tbls/cteq6/Makefile
cteq-tbls/cteq6m/Makefile
cteq-tbls/cteq6std/Makefile
cteq-tbls/cteq61/Makefile
cteq-tbls/ctq65-pds/Makefile
cteq-tbls/ctq65c/Makefile
cteq-tbls/ctq65s/Makefile
cteq-tbls/ctq66a/Makefile
cteq-tbls/ctq66c/Makefile
cteq-tbls/ctq66m/Makefile
])
AC_OUTPUT
