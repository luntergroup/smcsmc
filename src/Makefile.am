bin_PROGRAMS = smcsmc smcsmc_dbg
#smcsmc_omp smcsmc_dbg
#particleCountainer particle param coalevent scrm_dummy 
#check_PROGRAMS = smcsmc_dbg smcsmc_prof smcsmc_omp_prof unit_tests generate_vcf smcsmc_omp_prof smcsmc_omp
check_PROGRAMS = unit_tests  smcsmc_prof  generate_vcf

TESTS = unit_tests 

PROG = SMCSMC 

generate_vcf: 
	test -s ./sim-1Samples2msdata1.vcf  || ./unittest/pop_struct.py sim-1 2 1	
	test -s ./sim-1Samples2msdata1.rgvcf || ./unittest/gvcf2rgvcf.py -i sim-1Samples2msdata1.gvcf -seqlen 30000000
	test -s sim-1Samples4msdata1_edited.gvcf || { ./unittest/pop_struct.py sim-1 4 1 ; \
	head -50 sim-1Samples4msdata1.vcf > sim-1Samples4msdata1_first50_lines.vcf ; \
	sed -e "s/END=1434;/END=1430;/" -e "s/END=6582;/END=5000;/" -e "s/END=24916;/END=24900;/" -e "s/END=35798;/END=34000;/" sim-1Samples4msdata1.gvcf  | head -150 > sim-1Samples4msdata1_edited.gvcf ; \
	tail -n+3150 sim-1Samples4msdata1.gvcf >> sim-1Samples4msdata1_edited.gvcf; }
	test -s ./sim-1Samples4msdata1_edited.rgvcf || ./unittest/gvcf2rgvcf.py -i sim-1Samples4msdata1_edited.gvcf -seqlen 30000000

unittests: unit_tests 
#.PHONY: generate_vcf
	./unit_tests 

common_flags = -O3 -Wall -Wno-unknown-pragmas -I scrm/src/
smcsmc_CXXFLAGS= -std=c++0x -g $(common_flags) -DNDEBUG   #options for not showing debug messages
#smcsmc_CXXFLAGS= -std=c++0x -g $(common_flags) -DNDEBUG -fopenmp  #options for not showing debug messages

smcsmc_prof_CXXFLAGS= -std=c++0x -pg $(common_flags) -DNDEBUG   #options for not showing debug messages
smcsmc_dbg_CXXFLAGS= -std=c++0x -I scrm/src/ #options for showing debug messages

smcsmc_omp_CXXFLAGS= -std=c++0x -g $(common_flags) -DNDEBUG -fopenmp  #options for not showing debug messages
#smcsmc_omp_CXXFLAGS= -std=c++0x -g $(common_flags) -DNDEBUG   #options for not showing debug messages

smcsmc_omp_prof_CXXFLAGS= -std=c++0x -pg $(common_flags) -fopenmp -DNDEBUG   #options for not showing debug messages



#smcsmc_prof_CXXFLAGS= -std=c++0x -pg $(common_flags) -fopenmp -DNDEBUG   #options for not showing debug messages
#smcsmc_dbg_CXXFLAGS= -std=c++0x -fopenmp #options for showing debug messages
#scrm/src/summary_statistics/oriented_forest.cc 

scrm_param = scrm/src/param.cc scrm/src/param.h \
        scrm/src/model.cc scrm/src/model.h \
        scrm/src/summary_statistics/newick_tree.cc \
        scrm/src/summary_statistics/frequency_spectrum.cc \
        scrm/src/summary_statistics/seg_sites.cc \
		scrm/src/summary_statistics/tmrca.cc \
		scrm/src/summary_statistics/first_last_tmrca.cc 
		
		        

scrm_src = scrm/src/node.cc scrm/src/node.h \
        scrm/src/node_container.cc scrm/src/node_container.h \
        scrm/src/time_interval.cc scrm/src/time_interval.h \
        scrm/src/tree_point.cc scrm/src/tree_point.h \
        scrm/src/forest.cc \
        scrm/src/random/mersenne_twister.h  \
        scrm/src/random/random_generator.cc \
        scrm/src/random/mersenne_twister.cc   scrm/src/random/random_generator.h \
        scrm/src/random/fastfunc.cc           scrm/src/random/fastfunc.h 
        
#scrm_linear_SOURCES = $(scrm_param) $(scrm_src) scrm_omp.cpp
#scrm_linear_prof_SOURCES = $(scrm_param) $(scrm_src) scrm_omp.cpp
#scrm_omp_SOURCES = $(scrm_param) $(scrm_src) scrm_omp.cpp

debug_src = scrm/src/random/constant_generator.cc scrm/src/random/constant_generator.h \
            scrm/src/forest-debug.cc \
            debug/particle-debug.cpp \
            debug/count-debug.cpp 

arg_src = coalevent.cpp \
          particle.cpp \
          particleContainer.cpp 

count_src = count.cpp count-semi-debug.cpp 
       
param_src = pattern.cpp pfparam.cpp segdata.cpp rescue.cpp 

series_src = particleContainer_extra.cpp count_extra.cpp

omp_src = particleContainer_extra_omp.cpp count_extra_omp.cpp

shared_SOURCES = $(arg_src) $(scrm_param) $(scrm_src) $(param_src) $(count_src) help.cpp 

smcsmc_shared_SOURCES = $(shared_SOURCES) smcsmc.cpp

smcsmc_SOURCES = $(smcsmc_shared_SOURCES) debug/usage_ghost.cpp $(series_src) 
#smcsmc_SOURCES = $(shared_SOURCES) debug/usage_ghost.cpp $(omp_src) 

smcsmc_omp_SOURCES = $(smcsmc_shared_SOURCES)  debug/usage_ghost.cpp $(omp_src)
#smcsmc_omp_SOURCES = $(shared_SOURCES)  debug/usage_ghost.cpp count_extra_omp.cpp  particleContainer_extra.cpp
#smcsmc_omp_SOURCES = $(shared_SOURCES)  debug/usage_ghost.cpp count_extra.cpp  particleContainer_extra_omp.cpp

#smcsmc_prof_SOURCES = $(shared_SOURCES)  smcsmc.cpp debug/usage.cpp 
smcsmc_prof_SOURCES = $(smcsmc_SOURCES)  

smcsmc_omp_prof_SOURCES = $(smcsmc_omp_SOURCES)  


smcsmc_dbg_SOURCES = $(smcsmc_SOURCES) $(debug_src) 


unit_tests_SOURCES = $(shared_SOURCES) $(series_src)  \
                     $(debug_src) \
                     unittest/test_runner.cpp \
                     unittest/test_pattern.cpp \
                     unittest/test_Coalevent.cpp \
                     unittest/test_FroestState.cpp\
                     unittest/test_PfCount.cpp\
                     unittest/test_ParticleContainer.cpp
                                          
unit_tests_CXXFLAGS = -DUNITTEST -DNDEBUG -std=c++0x -I scrm/src/ -I src/

#unit_tests_LDADD= -L/home/joezhu/cppunit/cppunit-1.12.1/lib -lcppunit -ldl
unit_tests_LDADD= -lcppunit -ldl
# -L/opt/local/lib -lcppunit -ldl  #link the cppunit unittest library in mac, cppunit was installed via macports
#unit_tests_LDADD= -L/usr/include 

include_HEADERS = segdata.hpp pfparam.hpp particle.hpp pattern.hpp\
                  scrm/src/forest.h scrm/src/model.h scrm/src/node.h scrm/src/node_container.h \
                  scrm/src/time_interval.h scrm/src/tree_point.h \
                  scrm/src/random/constant_generator.h scrm/src/random/mersenne_twister.h \
                  scrm/src/random/random_generator.h scrm/src/event.h \
                  scrm/src/summary_statistics/newick_tree.h \
                  scrm/src/summary_statistics/frequency_spectrum.h \
                  scrm/src/summary_statistics/seg_sites.h \
                  scrm/src/summary_statistics/tmrca.h 
                  
#scrm/src/contemporaries_container.h scrm/src/summary_statistics/oriented_forest.h 
				  
clean-local: 
	test -z "unittest/sim-1Samples2msdata1.seg" || rm -f unittest/sim-1*
