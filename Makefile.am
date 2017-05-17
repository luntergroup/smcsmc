EXTRA_DIST = Doxyfile bootstrap pfpsmc/R/*r pfpsmc/R/ms_vcf.py #acinclude.m4
SMCSMCVERSION = $(shell git show HEAD | head -1 | sed -e "s/commit //g" | cat)
SCRMVERSION = $(shell git submodule status | grep scrm | sed -e "s/ //g" -e "s/src.*//" | cat)
COMPILEDATE = $(shell date -u | sed -e "s/ /-/g")

# disable for Darwin and rescomp
#TCMALLOC = -ltcmalloc

# set for rescomp
boostinclude = -I/apps/well/boost/1.59.0-gcc5.4.0-py27/include
boostldflags = -L/apps/well/boost/1.59.0-gcc5.4.0-py27/lib -Wl,-rpath,/apps/well/boost/1.59.0-gcc5.4.0-py27

#############################################

distdir = $(PACKAGE)-$(VERSION)

bin_PROGRAMS = smcsmc smcsmc_dbg smcsmc_static
check_PROGRAMS = unittests smcsmc_prof smcsmc_gperf
scrmcheck_PROGRAMS = scrm smcsmcSCRM smcsmcSCRM_dbg scrm_dbg
scrmcheckdir = bindir

PROG = SMCSMC

TESTS = unittests

.PHONY: scrmcheck
scrmcheck:
	make scrm scrm_dbg
	make smcsmcSCRM
	make smcsmcSCRM_dbg


.PHONY: generate_vcf$(EXEEXT)
generate_vcf$(EXEEXT):
	test -s ./sim-1Samples2msdata1.vcf  || ./unittest/pop_struct.py sim-1 2 1
	test -s ./sim-1Samples2msdata1.rgvcf || ./unittest/gvcf2rgvcf.py -i sim-1Samples2msdata1.gvcf -seqlen 30000000
	test -s sim-1Samples4msdata1_edited.gvcf || { ./unittest/pop_struct.py sim-1 4 1 ; \
	head -50 sim-1Samples4msdata1.vcf > sim-1Samples4msdata1_first50_lines.vcf ; \
	sed -e "s/END=1434;/END=1430;/" -e "s/END=6582;/END=5000;/" -e "s/END=24916;/END=24900;/" -e "s/END=35798;/END=34000;/" sim-1Samples4msdata1.gvcf  | head -150 > sim-1Samples4msdata1_edited.gvcf ; \
	tail -n+3150 sim-1Samples4msdata1.gvcf >> sim-1Samples4msdata1_edited.gvcf; }
	test -s ./sim-1Samples4msdata1_edited.rgvcf || ./unittest/gvcf2rgvcf.py -i sim-1Samples4msdata1_edited.gvcf -seqlen 30000000


common_flags = -std=c++11 -O3 -Wall -Wno-unknown-pragmas -Isrc/scrm/src/ -Isrc/scrm/src/random -DSMCSMCVERSION=\"${SMCSMCVERSION}\" -DSCRMVERSION=\"${SCRMVERSION}\" -DCOMPILEDATE=\"${COMPILEDATE}\" -DLANE_MAX_SIZE=50 -g

# Binary flags
smcsmc_CXXFLAGS          = $(common_flags) -DNDEBUG -O3 $(boostinclude)
smcsmcSCRM_CXXFLAGS      = $(smcsmc_CXXFLAGS) -D_SCRM
scrm_CXXFLAGS            = $(common_flags) -DNDEBUG -O3
smcsmc_LDADD             = -lboost_iostreams -lz ${TCMALLOC}
scrm_LDADD               = -lboost_iostreams -lz ${TCMALLOC}
smcsmc_LDFLAGS           = $(boostldflags)
scrm_LDFLAGS             = $(boostldflags)

# Debugger flags
scrm_dbg_CXXFLAGS            = $(common_flags)
smcsmc_dbg_CXXFLAGS          = -std=c++0x -I src/scrm/src/ -Isrc/ -Isrc/scrm/src/random -fno-omit-frame-pointer -O3
smcsmcSCRM_dbg_CXXFLAGS      = $(smcsmc_dbg_CXXFLAGS) -D_SCRM
smcsmc_dbg_LDADD             = -lboost_iostreams -lz ${TCMALLOC}

# Profiler flags
smcsmc_prof_CXXFLAGS     = $(common_flags) -DNDEBUG -fno-omit-frame-pointer -std=c++0x -pg -O3
smcsmc_prof_LDADD        = -lboost_iostreams -lz ${TCMALLOC}

# Static flags
smcsmc_static_CXXFLAGS   = $(smcsmc_CXXFLAGS)
smcsmc_static_LDFLAGS    = -static -static-libgcc -static-libstdc++
smcsmc_static_LDADD      = -lboost_iostreams -lz ${TCMALLOC}

# Gperf flags
smcsmc_gperf_CXXFLAGS     = $(common_flags) -DNDEBUG -fno-omit-frame-pointer -std=c++0x -DGPERF -O3
smcsmc_gperf_LDADD        = /usr/lib/libtcmalloc.so -lboost_iostreams -lz -lprofiler

# Notes:
#
# - run "sudo apt-get install google-perftools libgoogle-perftools-dev"
# - to compile, "make smcsmc_gperf"
# - to run, "env HEAPPROFILE=/full/path/to/profileoutput ./smcsmc_gperf [parameters]"
# - to analyze, "google-pprof --gv ./smcsmc_gperf /full/path/to/profileoutput.<number>.heap"

# Tester flags
unittests_CXXFLAGS       = -DUNITTEST -DNDEBUG -std=c++0x -Isrc/ -I src/scrm/src/ -Isrc/scrm/src/random --coverage
unittests_LDADD          = -lcppunit -ldl -lboost_iostreams -lz



scrm_src = src/scrm/src/node.cc \
           src/scrm/src/node.h \
           src/scrm/src/node_container.cc \
           src/scrm/src/node_container.h \
           src/scrm/src/time_interval.cc \
           src/scrm/src/time_interval.h \
           src/scrm/src/tree_point.cc src/scrm/src/tree_point.h \
           src/scrm/src/forest.cc \
           src/scrm/src/forest-debug.cc \
           src/scrm/src/random/mersenne_twister.h  \
           src/scrm/src/random/random_generator.cc \
           src/scrm/src/random/mersenne_twister.cc \
           src/scrm/src/random/random_generator.h \
           src/scrm/src/random/fastfunc.cc \
           src/scrm/src/random/fastfunc.h \
           src/scrm/src/param.cc \
           src/scrm/src/param.h \
           src/scrm/src/model.cc \
           src/scrm/src/model.h \
           src/scrm/src/summary_statistics/newick_tree.cc \
           src/scrm/src/summary_statistics/frequency_spectrum.cc \
           src/scrm/src/summary_statistics/seg_sites.cc \
           src/scrm/src/summary_statistics/tmrca.cc \
           src/scrm/src/summary_statistics/oriented_forest.cc

scrm_debug_src = src/scrm/src/random/constant_generator.cc src/scrm/src/random/constant_generator.h \
                 src/scrm/src/random/constant_generator.h

debug_src = src/scrm/src/random/constant_generator.cc \
            src/scrm/src/random/constant_generator.h

smcsmc_src = src/coalevent.cpp \
             src/particle.cpp \
             src/particleContainer.cpp \
             src/count.cpp \
             src/count-semi-debug.cpp \
             src/pattern.cpp \
             src/pfparam.cpp \
             src/segdata.cpp \
             src/arena.cpp \
             src/model_summary.cpp

shared_SOURCES = $(smcsmc_src) $(scrm_src)
smcsmc_shared_SOURCES = $(shared_SOURCES)

# Binary sources
smcsmc_SOURCES          = $(smcsmc_shared_SOURCES) src/smcsmc.cpp
smcsmcSCRM_SOURCES      = $(smcsmc_SOURCES)
scrm_SOURCES            = $(scrm_src) src/scrm/src/scrm.cc

# Debugger sources
scrm_dbg_SOURCES  = $(scrm_debug_src) $(scrm_SOURCES)
smcsmc_dbg_SOURCES = $(smcsmc_SOURCES) $(debug_src)
smcsmcSCRM_dbg_SOURCES = $(smcsmc_SOURCES) $(debug_src)

# Static sources
smcsmc_static_SOURCES = $(smcsmc_SOURCES)

# Profiler sources
smcsmc_prof_SOURCES = $(smcsmc_SOURCES)
smcsmc_gperf_SOURCES = $(smcsmc_SOURCES)

test_src = tests/unittest/test_runner.cpp \
		   tests/unittest/test_TwoDoubles.cpp \
		   tests/unittest/test_Coalevent.cpp \
		   tests/unittest/test_FroestState.cpp \
		   tests/unittest/test_ParticleContainer.cpp \
		   tests/unittest/test_Pattern.cpp \
		   tests/unittest/test_PfParam.cpp \
		   tests/unittest/test_Segment.cpp

unittests_SOURCES = $(smcsmc_shared_SOURCES)  $(test_src)

.PHONY: clean-local
clean-local:
	test -z "unittest/sim-1Samples2msdata1.seg" || rm -f unittest/sim-1*
