bin_PROGRAMS = pf-ARG pf-ARG_omp pf-ARG_dbg
#particleCountainer particle param coalevent scrm_dummy
check_PROGRAMS = pf-ARG_dbg pf-ARG_prof pf-ARG_omp_prof unit_tests 
TESTS = unit_tests 

PROG = PF-ARG 

pf_ARG_CXXFLAGS= -std=c++0x -g -O3 -Wall -Wno-unknown-pragmas -DNDEBUG   #options for not showing debug messages
#pf_ARG_CXXFLAGS= -std=c++0x -g -O3 -Wall -Wno-unknown-pragmas -DNDEBUG -fopenmp  #options for not showing debug messages

pf_ARG_prof_CXXFLAGS= -std=c++0x -pg -O3 -Wall -Wno-unknown-pragmas  -DNDEBUG   #options for not showing debug messages
pf_ARG_dbg_CXXFLAGS= -std=c++0x  #options for showing debug messages

pf_ARG_omp_CXXFLAGS= -std=c++0x -g -O3 -Wall -Wno-unknown-pragmas -DNDEBUG -fopenmp  #options for not showing debug messages
#pf_ARG_omp_CXXFLAGS= -std=c++0x -g -O3 -Wall -Wno-unknown-pragmas -DNDEBUG   #options for not showing debug messages

pf_ARG_omp_prof_CXXFLAGS= -std=c++0x -pg -O3 -Wall -Wno-unknown-pragmas -fopenmp -DNDEBUG   #options for not showing debug messages



#pf_ARG_prof_CXXFLAGS= -std=c++0x -pg -O3 -Wall -Wno-unknown-pragmas -fopenmp -DNDEBUG   #options for not showing debug messages
#pf_ARG_dbg_CXXFLAGS= -std=c++0x -fopenmp #options for showing debug messages



scrm_param = scrm/param.cc scrm/param.h \
        scrm/model.cc scrm/model.h \
        scrm/summary_statistics/newick_tree.cc \
        scrm/summary_statistics/frequency_spectrum.cc \
        scrm/summary_statistics/seg_sites.cc \
		scrm/summary_statistics/tmrca.cc \
		scrm/summary_statistics/first_last_tmrca.cc
		        

scrm_src = scrm/node.cc scrm/node.h \
        scrm/node_container.cc scrm/node_container.h \
        scrm/time_interval.cc scrm/time_interval.h \
        scrm/tree_point.cc scrm/tree_point.h \
        scrm/forest.cc \
        scrm/random/mersenne_twister.h  \
        scrm/random/random_generator.cc \
        scrm/random/mersenne_twister.cc   scrm/random/random_generator.h \
        scrm/random/fastfunc.cc           scrm/random/fastfunc.h 
        
#scrm_linear_SOURCES = $(scrm_param) $(scrm_src) scrm_omp.cpp
#scrm_linear_prof_SOURCES = $(scrm_param) $(scrm_src) scrm_omp.cpp
#scrm_omp_SOURCES = $(scrm_param) $(scrm_src) scrm_omp.cpp

debug_src = scrm/random/constant_generator.cc scrm/random/constant_generator.h \
            scrm/forest-debug.cc \
            debug/particle-debug.cpp \
            debug/count-debug.cpp \
            debug/variantReader-debug.cpp


arg_src = coalevent.cpp \
          particle.cpp \
          particleContainer.cpp 

count_src = count.cpp count-semi-debug.cpp 
       
param_src = pattern.cpp pfparam.cpp variantReader.cpp rescue.cpp

series_src = particleContainer_extra.cpp count_extra.cpp

omp_src = particleContainer_extra_omp.cpp count_extra_omp.cpp

shared_SOURCES = $(arg_src) $(scrm_param) $(scrm_src) $(param_src) $(count_src) pf-ARG.cpp

pf_ARG_SOURCES = $(shared_SOURCES) debug/usage_ghost.cpp $(series_src) 
#pf_ARG_SOURCES = $(shared_SOURCES) debug/usage_ghost.cpp $(omp_src) 

pf_ARG_omp_SOURCES = $(shared_SOURCES)  debug/usage_ghost.cpp $(omp_src)
#pf_ARG_omp_SOURCES = $(shared_SOURCES)  debug/usage_ghost.cpp count_extra_omp.cpp  particleContainer_extra.cpp
#pf_ARG_omp_SOURCES = $(shared_SOURCES)  debug/usage_ghost.cpp count_extra.cpp  particleContainer_extra_omp.cpp

#pf_ARG_prof_SOURCES = $(shared_SOURCES)  pf-ARG.cpp debug/usage.cpp 
pf_ARG_prof_SOURCES = $(pf_ARG_SOURCES)  

pf_ARG_omp_prof_SOURCES = $(pf_ARG_omp_SOURCES)  


pf_ARG_dbg_SOURCES = $(pf_ARG_SOURCES) $(debug_src) 


unit_tests_SOURCES = $(arg_src) $(scrm_param) $(scrm_src) $(param_src) $(count_src) $(series_src)  \
                     $(debug_src) \
                     tests/test_runner.cpp \
                     tests/test_param.cpp \
                     tests/test_pattern.cpp \
                     tests/test_Coalevent.cpp \
                     tests/test_FroestState.cpp\
                     tests/test_PfCount.cpp\
                     tests/test_ParticleContainer.cpp\
                     tests/test_variantReader.cpp \
                     #tests/test_random_generator.cpp \
		     #tests/test_Count.cpp\
                     #$(pf_ARG_SOURCES)

unit_tests_CXXFLAGS = -DUNITTEST -DNDEBUG -std=c++0x 

unit_tests_LDADD= -L/home/joezhu/cppunit/cppunit-1.12.1/lib -lcppunit -ldl
# -L/opt/local/lib -lcppunit -ldl  #link the cppunit unittest library in mac, cppunit was installed via macports
#unit_tests_LDADD= -L/usr/include 

include_HEADERS = variantReader.hpp pfparam.hpp particle.hpp pattern.hpp\
                  scrm/forest.h scrm/model.h scrm/node.h scrm/node_container.h \
                  scrm/time_interval.h scrm/tree_point.h \
                  scrm/random/constant_generator.h scrm/random/mersenne_twister.h \
                  scrm/random/random_generator.h scrm/event.h \
                  scrm/summary_statistics/newick_tree.h \
                  scrm/summary_statistics/frequency_spectrum.h \
                  scrm/summary_statistics/seg_sites.h \
                  scrm/summary_statistics/tmrca.h
				  
