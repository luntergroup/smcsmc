EXTRA_DIST = Doxyfile bootstrap pfpsmc/R/*r pfpsmc/R/ms_vcf.py #acinclude.m4
distdir = $(PACKAGE)-$(VERSION)
#CXXFLAGS= -I/opt/local/include
#LDFLAGS= -L/opt/local/lib # this is used to link to the unittest library

#include doxygen-include.am
#include aminclude.am 

bin_PROGRAMS = smcsmc smcsmc_dbg smcsmcRECOMBOFF smcsmcSCRM smcsmcRECOMBOFF_dbg
#smcsmc_omp smcsmc_dbg smcsmcRECOMBOFF_dbg  smcsmcSCRM
#particleCountainer particle param coalevent scrm_dummy 
#check_PROGRAMS = smcsmc_dbg smcsmc_prof smcsmc_omp_prof unittests generate_vcf smcsmc_omp_prof smcsmc_omp
check_PROGRAMS = unittests  smcsmc_prof scrm

TESTS = unittests 

PROG = SMCSMC 

.PHONY: generate_vcf$(EXEEXT) 
generate_vcf$(EXEEXT): 
	test -s ./sim-1Samples2msdata1.vcf  || ./unittest/pop_struct.py sim-1 2 1	
	test -s ./sim-1Samples2msdata1.rgvcf || ./unittest/gvcf2rgvcf.py -i sim-1Samples2msdata1.gvcf -seqlen 30000000
	test -s sim-1Samples4msdata1_edited.gvcf || { ./unittest/pop_struct.py sim-1 4 1 ; \
	head -50 sim-1Samples4msdata1.vcf > sim-1Samples4msdata1_first50_lines.vcf ; \
	sed -e "s/END=1434;/END=1430;/" -e "s/END=6582;/END=5000;/" -e "s/END=24916;/END=24900;/" -e "s/END=35798;/END=34000;/" sim-1Samples4msdata1.gvcf  | head -150 > sim-1Samples4msdata1_edited.gvcf ; \
	tail -n+3150 sim-1Samples4msdata1.gvcf >> sim-1Samples4msdata1_edited.gvcf; }
	test -s ./sim-1Samples4msdata1_edited.rgvcf || ./unittest/gvcf2rgvcf.py -i sim-1Samples4msdata1_edited.gvcf -seqlen 30000000


##.PHONY: exec_unittests
##exec_unittests: unittests 
##./unittests 


common_flags = -std=c++0x -O3 -Wall -Wno-unknown-pragmas -Isrc/scrm/src/

# -DNDEBUG  # Define macro NDEBUG, for not showing debug messages
# -D_RecombRecordingOff # Define macro D_RecombRecordingOff, for turning of Recombevent recorders and recombination rate update

# Binary flags
smcsmc_CXXFLAGS          = $(common_flags) -DNDEBUG -O3
smcsmcRECOMBOFF_CXXFLAGS = $(smcsmc_CXXFLAGS) -D_RecombRecordingOff
smcsmcSCRM_CXXFLAGS      = $(smcsmcRECOMBOFF_CXXFLAGS) -D_SCRM
smcsmc_omp_CXXFLAGS      = $(common_flags) -DNDEBUG -fopenmp -O3
scrm_CXXFLAGS            = $(common_flags) -DNDEBUG -O3

# Debugger flags
smcsmc_dbg_CXXFLAGS          = -std=c++0x -I src/scrm/src/ -Isrc/ -fno-omit-frame-pointer -O2
smcsmcRECOMBOFF_dbg_CXXFLAGS = $(smcsmc_dbg_CXXFLAGS)

# Profiler flags
smcsmc_prof_CXXFLAGS     = $(common_flags) -DNDEBUG -fno-omit-frame-pointer -std=c++0x -pg -O1
smcsmc_omp_prof_CXXFLAGS = $(smcsmc_prof_CXXFLAGS) -fopenmp 

# Tester flags
unittests_CXXFLAGS       = -DUNITTEST -DNDEBUG -std=c++0x -Isrc/ -I src/scrm/src/
unittests_LDADD          = -lcppunit -ldl

##export SMCSMCROOT=src
##export SCRMROOT=src/scrm/src
##export UNITTESTROOT=tests/unittest

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
           src/scrm/src/summary_statistics/first_last_tmrca.cc

#scrm_linear_SOURCES = $(scrm_param) $(scrm_src) scrm_omp.cpp
#scrm_linear_prof_SOURCES = $(scrm_param) $(scrm_src) scrm_omp.cpp
#scrm_omp_SOURCES = $(scrm_param) $(scrm_src) scrm_omp.cpp

debug_src = src/scrm/src/random/constant_generator.cc \
            src/scrm/src/random/constant_generator.h \
            src/debug/particle-debug.cpp \
            src/debug/count-debug.cpp 

smcsmc_src = src/coalevent.cpp \
             src/particle.cpp \
             src/particleContainer.cpp \
             src/count.cpp \
             src/count-semi-debug.cpp \
             src/pattern.cpp \
             src/pfparam.cpp \
             src/segdata.cpp \
             src/rescue.cpp \
             src/help.cpp \
             src/arena.cpp

series_src = src/particleContainer_extra.cpp \
             src/count_extra.cpp

omp_src = src/particleContainer_extra_omp.cpp \
          src/count_extra_omp.cpp


shared_SOURCES = $(smcsmc_src) $(scrm_src)
smcsmc_shared_SOURCES = $(shared_SOURCES) 

# Binary sources
smcsmc_SOURCES          = $(smcsmc_shared_SOURCES) src/debug/usage_ghost.cpp $(series_src) src/smcsmc.cpp
smcsmcRECOMBOFF_SOURCES = $(smcsmc_SOURCES)
smcsmcSCRM_SOURCES      = $(smcsmc_SOURCES)
smcsmc_omp_SOURCES      = $(smcsmc_shared_SOURCES) src/debug/usage_ghost.cpp $(omp_src)
scrm_SOURCES            = $(scrm_src) src/scrm/src/scrm.cc

# Debugger sources
smcsmc_dbg_SOURCES = $(smcsmc_SOURCES) $(debug_src) 
smcsmcRECOMBOFF_dbg_SOURCES = $(smcsmc_dbg_SOURCES)

# Profiler sources
smcsmc_prof_SOURCES = $(smcsmc_SOURCES)  
smcsmc_omp_prof_SOURCES = $(smcsmc_omp_SOURCES)  

#smcsmc_omp_SOURCES = $(shared_SOURCES)  debug/usage_ghost.cpp count_extra_omp.cpp  particleContainer_extra.cpp
#smcsmc_omp_SOURCES = $(shared_SOURCES)  debug/usage_ghost.cpp count_extra.cpp  particleContainer_extra_omp.cpp
#smcsmc_prof_SOURCES = $(shared_SOURCES)  smcsmc.cpp debug/usage.cpp 

# Tester sources
test_src = tests/unittest/test_pattern.cpp \
	   tests/unittest/test_runner.cpp \
	   tests/unittest/test_Coalevent.cpp \
           tests/unittest/test_FroestState.cpp\
           tests/unittest/test_PfCount.cpp \
           tests/unittest/test_ParticleContainer.cpp

unittests_SOURCES = $(smcsmc_shared_SOURCES)  $(test_src) $(series_src)

# -L/opt/local/lib -lcppunit -ldl  #link the cppunit unittest library in mac, cppunit was installed via macports
#-Isrc/scrm/src/ -Isrc/ -Itests/unittest/
 #including hearders for make distribution
#include_HEADERS = segdata.hpp pfparam.hpp particle.hpp pattern.hpp\
                  #src/scrm/src/forest.h src/scrm/src/model.h src/scrm/src/node.h src/scrm/src/node_container.h \
                  #src/scrm/src/time_interval.h src/scrm/src/tree_point.h \
                  #src/scrm/src/random/constant_generator.h src/scrm/src/random/mersenne_twister.h \
                  #src/scrm/src/random/random_generator.h src/scrm/src/event.h \
                  #src/scrm/src/summary_statistics/newick_tree.h \
                  #src/scrm/src/summary_statistics/frequency_spectrum.h \
                  #src/scrm/src/summary_statistics/seg_sites.h \
                  #src/scrm/src/summary_statistics/tmrca.h 

.PHONY: clean-local                  
clean-local: 
	test -z "unittest/sim-1Samples2msdata1.seg" || rm -f unittest/sim-1*

