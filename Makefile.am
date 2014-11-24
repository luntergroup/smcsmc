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

common_flags = -O3 -Wall -Wno-unknown-pragmas -I $(scrm_root)/
smcsmc_CXXFLAGS= -std=c++0x -g $(common_flags) -DNDEBUG   


smcsmc_prof_CXXFLAGS= -std=c++0x -pg $(common_flags) -DNDEBUG   #options for not showing debug messages
smcsmc_dbg_CXXFLAGS= -std=c++0x -I $(scrm_root)/ -I$(smcsmc_root)/ #options for showing debug messages

smcsmc_omp_CXXFLAGS= -std=c++0x -g $(common_flags) -DNDEBUG -fopenmp  #options for not showing debug messages
#smcsmc_omp_CXXFLAGS= -std=c++0x -g $(common_flags) -DNDEBUG   #options for not showing debug messages

smcsmc_omp_prof_CXXFLAGS= -std=c++0x -pg $(common_flags) -fopenmp -DNDEBUG   #options for not showing debug messages



#smcsmc_prof_CXXFLAGS= -std=c++0x -pg $(common_flags) -fopenmp -DNDEBUG   #options for not showing debug messages
#smcsmc_dbg_CXXFLAGS= -std=c++0x -fopenmp #options for showing debug messages


smcsmc_root = src

scrm_root = src/scrm/src

scrm_src = $(scrm_root)/node.cc \
           $(scrm_root)/node.h \
           $(scrm_root)/node_container.cc \
           $(scrm_root)/node_container.h \
           $(scrm_root)/time_interval.cc \
           $(scrm_root)/time_interval.h \
           $(scrm_root)/tree_point.cc $(scrm_root)/tree_point.h \
           $(scrm_root)/forest.cc \
           $(scrm_root)/random/mersenne_twister.h  \
           $(scrm_root)/random/random_generator.cc \
           $(scrm_root)/random/mersenne_twister.cc \
           $(scrm_root)/random/random_generator.h \
           $(scrm_root)/random/fastfunc.cc \
           $(scrm_root)/random/fastfunc.h \
           $(scrm_root)/param.cc \
           $(scrm_root)/param.h \
           $(scrm_root)/model.cc \
           $(scrm_root)/model.h \
           $(scrm_root)/summary_statistics/newick_tree.cc \
           $(scrm_root)/summary_statistics/frequency_spectrum.cc \
           $(scrm_root)/summary_statistics/seg_sites.cc \
		   $(scrm_root)/summary_statistics/tmrca.cc \
		   $(scrm_root)/summary_statistics/first_last_tmrca.cc 
        
#scrm_linear_SOURCES = $(scrm_param) $(scrm_src) scrm_omp.cpp
#scrm_linear_prof_SOURCES = $(scrm_param) $(scrm_src) scrm_omp.cpp
#scrm_omp_SOURCES = $(scrm_param) $(scrm_src) scrm_omp.cpp

debug_src = $(scrm_root)/random/constant_generator.cc \
            $(scrm_root)/random/constant_generator.h \
            $(scrm_root)/forest-debug.cc \
            $(smcsmc_root)/debug/particle-debug.cpp \
            $(smcsmc_root)/debug/count-debug.cpp 

smcsmc_src = $(smcsmc_root)/coalevent.cpp \
             $(smcsmc_root)/particle.cpp \
             $(smcsmc_root)/particleContainer.cpp \
             $(smcsmc_root)/count.cpp \
             $(smcsmc_root)/count-semi-debug.cpp \
             $(smcsmc_root)/pattern.cpp \
             $(smcsmc_root)/pfparam.cpp \
             $(smcsmc_root)/segdata.cpp \
             $(smcsmc_root)/rescue.cpp \
             $(smcsmc_root)/help.cpp

series_src = $(smcsmc_root)/particleContainer_extra.cpp \
             $(smcsmc_root)/count_extra.cpp

omp_src = $(smcsmc_root)/particleContainer_extra_omp.cpp \
          $(smcsmc_root)/count_extra_omp.cpp

shared_SOURCES = $(smcsmc_src) $(scrm_src)

smcsmc_shared_SOURCES = $(shared_SOURCES) $(smcsmc_root)/smcsmc.cpp

smcsmc_SOURCES = $(smcsmc_shared_SOURCES) $(smcsmc_root)/debug/usage_ghost.cpp $(series_src) 
#smcsmc_SOURCES = $(shared_SOURCES) debug/usage_ghost.cpp $(omp_src) 

smcsmc_omp_SOURCES = $(smcsmc_shared_SOURCES) $(smcsmc_root)/debug/usage_ghost.cpp $(omp_src)
#smcsmc_omp_SOURCES = $(shared_SOURCES)  debug/usage_ghost.cpp count_extra_omp.cpp  particleContainer_extra.cpp
#smcsmc_omp_SOURCES = $(shared_SOURCES)  debug/usage_ghost.cpp count_extra.cpp  particleContainer_extra_omp.cpp

#smcsmc_prof_SOURCES = $(shared_SOURCES)  smcsmc.cpp debug/usage.cpp 
smcsmc_prof_SOURCES = $(smcsmc_SOURCES)  

smcsmc_omp_prof_SOURCES = $(smcsmc_omp_SOURCES)  

smcsmc_dbg_SOURCES = $(smcsmc_SOURCES) $(debug_src) 

unittest_root = tests/unittest

unittests_SOURCES = $(smcsmc_root)/pattern.cpp \
					$(smcsmc_root)/pattern.hpp \
                    $(unittest_root)/test_pattern.cpp \
                    $(scrm_root)/model.cc \
					$(unittest_root)/test_runner.cpp \
					$(unittest_root)/test_Coalevent.cpp \
					$(smcsmc_root)/particle.cpp \
					$(smcsmc_root)/particle.hpp \
                    $(unittest_root)/test_FroestState.cpp\
                    $(unittest_root)/test_PfCount.cpp \
                    $(unittest_root)/test_ParticleContainer.cpp
			  

unittests_CXXFLAGS = -DUNITTEST -DNDEBUG -std=c++0x -I$(smcsmc_root)/ -I $(scrm_root)/

unittests_LDADD= -lcppunit -ldl
# -L/opt/local/lib -lcppunit -ldl  #link the cppunit unittest library in mac, cppunit was installed via macports
#-I$(scrm_root)/ -I$(smcsmc_root)/ -Itests/unittest/
 #including hearders for make distribution
#include_HEADERS = segdata.hpp pfparam.hpp particle.hpp pattern.hpp\
                  #$(scrm_root)/forest.h $(scrm_root)/model.h $(scrm_root)/node.h $(scrm_root)/node_container.h \
                  #$(scrm_root)/time_interval.h $(scrm_root)/tree_point.h \
                  #$(scrm_root)/random/constant_generator.h $(scrm_root)/random/mersenne_twister.h \
                  #$(scrm_root)/random/random_generator.h $(scrm_root)/event.h \
                  #$(scrm_root)/summary_statistics/newick_tree.h \
                  #$(scrm_root)/summary_statistics/frequency_spectrum.h \
                  #$(scrm_root)/summary_statistics/seg_sites.h \
                  #$(scrm_root)/summary_statistics/tmrca.h 

.PHONY: clean-local                  
clean-local: 
	test -z "unittest/sim-1Samples2msdata1.seg" || rm -f unittest/sim-1*

