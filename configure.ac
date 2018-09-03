AC_INIT([smcsmc],[1.0],[gerton.lunter@well.ox.ac.uk])

# remove default -O2; we will set this by hand (or not) later
# also include -fno-omit-frame-pointer for better debugging/profiling
CXXFLAGS="-g -fno-omit-frame-pointer"

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
AC_CHECK_LIB(cppunit,TestCase,[])
# Checks for header files for scrm.
AC_HEADER_STDC
AC_LANG(C++) 
AC_CHECK_HEADERS([cassert cstddef cfloat stdexcept iostream vector map iomanip fstream ctime cmath boost/assign/std/vector.hpp boost/random.hpp])
AC_CHECK_HEADERS([boost/random.hpp boost/assign/std/vector.hpp],[],AC_MSG_WARN([User define 'boost' location. e.g. CXXFLAGS=-I/opt/local/include]))

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for libraries.
AC_CHECK_LIB([m], [log])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

AC_CHECK_FUNCS([pow])
AC_HEADER_STDBOOL

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
