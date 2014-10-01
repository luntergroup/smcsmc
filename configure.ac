AC_INIT([pf-ARG],esyscmd([sh -c "date \"+%Y%m%d\" | tr -d '\n'"]),[joezhu@well.ox.ac.uk])

#AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE([subdir-objects])

# Suppress output at compilation
m4_ifndef([AM_SILENT_RULES], [m4_define([AM_SILENT_RULES],[])])
AM_SILENT_RULES([yes])

# Load macros in 'm4'-dir
AC_CONFIG_MACRO_DIR([m4])

AC_PROG_INSTALL
AC_PREREQ
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_RANLIB

# Check for C++11
AX_CXX_COMPILE_STDCXX_11([noext])

#Checks for libraries
AC_CHECK_LIB(cppunit,TestCase,[],LDFLAGS=-L/opt/local/lib)
#AC_CHECK_LIB(cppunit,TestCase,[],LDFLAGS=-L/home/joezhu/cppunit/lib)
# Checks for header files for scrm.
AC_HEADER_STDC
AC_LANG(C++) 
#AC_CHECK_HEADERS([cassert cstddef cfloat stdexcept iostream vector map iomanip fstream ctime cmath boost/assign/std/vector.hpp boost/random.hpp])
#AC_CHECK_HEADERS([boost/random.hpp boost/assign/std/vector.hpp],[],AC_MSG_WARN([User define 'boost' location. e.g. CXXFLAGS=-I/opt/local/include]))
#AC_CHECK_HEADERS([boost/random.hpp boost/assign/std/vector.hpp],[],CXXFLAGS=-I../opt/local/include)


# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T


# Checks for libraries.
#AC_CHECK_LIB([m], [log])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

AC_CHECK_FUNCS([pow])
#AC_HEADER_STDBOOL

m4_include([m4/ac_doxygen.m4])
DX_INIT_DOXYGEN(pf-ARG, Doxyfile, html/)
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)


# Check for presence of pdfLaTeX
AC_CHECK_PROG(PDFLATEX, pdflatex, pdflatex)
AC_CHECK_PROG(BIBTEX, bibtex,bibtex)
if test -z "$PDFLATEX"; then
  AC_MSG_WARN([Unable to create PDF version of the user manual.])
fi

AM_CONDITIONAL([HAVE_PDFLATEX], test -n "$PDFLATEX")

AC_CONFIG_FILES([Makefile src/Makefile ])
AC_OUTPUT
#
