EXTRA_DIST = Doxyfile bootstrap pfpsmc/R/*r pfpsmc/R/ms_vcf.py #acinclude.m4
distdir = $(PACKAGE)-$(VERSION)
#CXXFLAGS= -I/opt/local/include
#LDFLAGS= -L/opt/local/lib # this is used to link to the unittest library

#include doxygen-include.am
#include aminclude.am 


bin_PROGRAMS = smcsmc smcsmc_dbg
#smcsmc_omp smcsmc_dbg
#particleCountainer particle param coalevent scrm_dummy 
#check_PROGRAMS = smcsmc_dbg smcsmc_prof smcsmc_omp_prof unittests generate_vcf smcsmc_omp_prof smcsmc_omp
check_PROGRAMS = unittests  smcsmc_prof 

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

## -DNDEBUG  options for not showing debug messages

common_flags = -O3 -Wall -Wno-unknown-pragmas -I $(SCRMROOT)/
smcsmc_CXXFLAGS= -std=c++0x -g $(common_flags) -DNDEBUG   


smcsmc_prof_CXXFLAGS= -std=c++0x -pg $(common_flags) -DNDEBUG   #options for not showing debug messages
smcsmc_dbg_CXXFLAGS= -std=c++0x -I $(SCRMROOT)/ -I$(SMCSMCROOT)/ #options for showing debug messages

smcsmc_omp_CXXFLAGS= -std=c++0x -g $(common_flags) -DNDEBUG -fopenmp  #options for not showing debug messages
#smcsmc_omp_CXXFLAGS= -std=c++0x -g $(common_flags) -DNDEBUG   #options for not showing debug messages

smcsmc_omp_prof_CXXFLAGS= -std=c++0x -pg $(common_flags) -fopenmp -DNDEBUG   #options for not showing debug messages



#smcsmc_prof_CXXFLAGS= -std=c++0x -pg $(common_flags) -fopenmp -DNDEBUG   #options for not showing debug messages
#smcsmc_dbg_CXXFLAGS= -std=c++0x -fopenmp #options for showing debug messages


SMCSMCROOT = src

SCRMROOT = src/scrm/src

scrm_src = $(SCRMROOT)/node.cc \
           $(SCRMROOT)/node.h \
           $(SCRMROOT)/node_container.cc \
           $(SCRMROOT)/node_container.h \
           $(SCRMROOT)/time_interval.cc \
           $(SCRMROOT)/time_interval.h \
           $(SCRMROOT)/tree_point.cc $(SCRMROOT)/tree_point.h \
           $(SCRMROOT)/forest.cc \
           $(SCRMROOT)/random/mersenne_twister.h  \
           $(SCRMROOT)/random/random_generator.cc \
           $(SCRMROOT)/random/mersenne_twister.cc \
           $(SCRMROOT)/random/random_generator.h \
           $(SCRMROOT)/random/fastfunc.cc \
           $(SCRMROOT)/random/fastfunc.h \
           $(SCRMROOT)/param.cc \
           $(SCRMROOT)/param.h \
           $(SCRMROOT)/model.cc \
           $(SCRMROOT)/model.h \
           $(SCRMROOT)/summary_statistics/newick_tree.cc \
           $(SCRMROOT)/summary_statistics/frequency_spectrum.cc \
           $(SCRMROOT)/summary_statistics/seg_sites.cc \
		   $(SCRMROOT)/summary_statistics/tmrca.cc \
		   $(SCRMROOT)/summary_statistics/first_last_tmrca.cc 
        
#scrm_linear_SOURCES = $(scrm_param) $(scrm_src) scrm_omp.cpp
#scrm_linear_prof_SOURCES = $(scrm_param) $(scrm_src) scrm_omp.cpp
#scrm_omp_SOURCES = $(scrm_param) $(scrm_src) scrm_omp.cpp

debug_src = $(SCRMROOT)/random/constant_generator.cc \
            $(SCRMROOT)/random/constant_generator.h \
            $(SCRMROOT)/forest-debug.cc \
            $(SMCSMCROOT)/debug/particle-debug.cpp \
            $(SMCSMCROOT)/debug/count-debug.cpp 

smcsmc_src = $(SMCSMCROOT)/coalevent.cpp \
             $(SMCSMCROOT)/particle.cpp \
             $(SMCSMCROOT)/particleContainer.cpp \
             $(SMCSMCROOT)/count.cpp \
             $(SMCSMCROOT)/count-semi-debug.cpp \
             $(SMCSMCROOT)/pattern.cpp \
             $(SMCSMCROOT)/pfparam.cpp \
             $(SMCSMCROOT)/segdata.cpp \
             $(SMCSMCROOT)/rescue.cpp \
             $(SMCSMCROOT)/help.cpp

series_src = $(SMCSMCROOT)/particleContainer_extra.cpp \
             $(SMCSMCROOT)/count_extra.cpp

omp_src = $(SMCSMCROOT)/particleContainer_extra_omp.cpp \
          $(SMCSMCROOT)/count_extra_omp.cpp

shared_SOURCES = $(smcsmc_src) $(scrm_src)

smcsmc_shared_SOURCES = $(shared_SOURCES) $(SMCSMCROOT)/smcsmc.cpp

smcsmc_SOURCES = $(smcsmc_shared_SOURCES) $(SMCSMCROOT)/debug/usage_ghost.cpp $(series_src) 
#smcsmc_SOURCES = $(shared_SOURCES) debug/usage_ghost.cpp $(omp_src) 

smcsmc_omp_SOURCES = $(smcsmc_shared_SOURCES) $(SMCSMCROOT)/debug/usage_ghost.cpp $(omp_src)
#smcsmc_omp_SOURCES = $(shared_SOURCES)  debug/usage_ghost.cpp count_extra_omp.cpp  particleContainer_extra.cpp
#smcsmc_omp_SOURCES = $(shared_SOURCES)  debug/usage_ghost.cpp count_extra.cpp  particleContainer_extra_omp.cpp

#smcsmc_prof_SOURCES = $(shared_SOURCES)  smcsmc.cpp debug/usage.cpp 
smcsmc_prof_SOURCES = $(smcsmc_SOURCES)  

smcsmc_omp_prof_SOURCES = $(smcsmc_omp_SOURCES)  

smcsmc_dbg_SOURCES = $(smcsmc_SOURCES) $(debug_src) 

UNITTESTROOT = tests/unittest

test_src = $(UNITTESTROOT)/test_pattern.cpp \
		   $(UNITTESTROOT)/test_runner.cpp \
		   $(UNITTESTROOT)/test_Coalevent.cpp \
           $(UNITTESTROOT)/test_FroestState.cpp\
           $(UNITTESTROOT)/test_PfCount.cpp \
           $(UNITTESTROOT)/test_ParticleContainer.cpp
		
unittests_SOURCES = $(SMCSMCROOT)/pattern.cpp \
					$(SMCSMCROOT)/pattern.hpp \      
                    $(SCRMROOT)/model.cc \
					$(SMCSMCROOT)/particle.cpp \
					$(SMCSMCROOT)/particle.hpp \
					$(test_src)

unittests_CXXFLAGS = -DUNITTEST -DNDEBUG -std=c++0x -I$(SMCSMCROOT)/ -I $(SCRMROOT)/

unittests_LDADD= -lcppunit -ldl
# -L/opt/local/lib -lcppunit -ldl  #link the cppunit unittest library in mac, cppunit was installed via macports
#-I$(SCRMROOT)/ -I$(SMCSMCROOT)/ -Itests/unittest/
 #including hearders for make distribution
#include_HEADERS = segdata.hpp pfparam.hpp particle.hpp pattern.hpp\
                  #$(SCRMROOT)/forest.h $(SCRMROOT)/model.h $(SCRMROOT)/node.h $(SCRMROOT)/node_container.h \
                  #$(SCRMROOT)/time_interval.h $(SCRMROOT)/tree_point.h \
                  #$(SCRMROOT)/random/constant_generator.h $(SCRMROOT)/random/mersenne_twister.h \
                  #$(SCRMROOT)/random/random_generator.h $(SCRMROOT)/event.h \
                  #$(SCRMROOT)/summary_statistics/newick_tree.h \
                  #$(SCRMROOT)/summary_statistics/frequency_spectrum.h \
                  #$(SCRMROOT)/summary_statistics/seg_sites.h \
                  #$(SCRMROOT)/summary_statistics/tmrca.h 

.PHONY: clean-local                  
clean-local: 
	test -z "unittest/sim-1Samples2msdata1.seg" || rm -f unittest/sim-1*

