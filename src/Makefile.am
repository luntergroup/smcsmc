bin_PROGRAMS = pf-ARG pf-ARG_dbg 
#particleCountainer particle param coalevent scrm_dummy
check_PROGRAMS = pf-ARG_dbg pf-ARG_prof unit_tests
TESTS = unit_tests 

PROG = PF-ARG 

pf_ARG_CXXFLAGS= -std=c++0x -g -O3 -Wall -Wno-unknown-pragmas -DNDEBUG   #options for not showing debug messages
pf_ARG_prof_CXXFLAGS= -std=c++0x -pg -O3 -Wall -Wno-unknown-pragmas -DNDEBUG   #options for not showing debug messages
pf_ARG_dbg_CXXFLAGS= -std=c++0x  #options for showing debug messages

#particleCountainer_CXXFLAGS = -std=c++0x -pg -O3 -Wall -Wno-unknown-pragmas -DNDEBUG   #options for not showing debug messages
#particle_CXXFLAGS = -std=c++0x -pg -O3 -Wall -Wno-unknown-pragmas -DNDEBUG   #options for not showing debug messages
#param_CXXFLAGS = -std=c++0x -pg -O3 -Wall -Wno-unknown-pragmas -DNDEBUG   #options for not showing debug messages
#coalevent_CXXFLAGS = -std=c++0x -pg -O3 -Wall -Wno-unknown-pragmas -DNDEBUG   #options for not showing debug messages
#scrm_dummy_CXXFLAGS = -std=c++0x -pg -O3 -Wall -Wno-unknown-pragmas -DNDEBUG   #options for not showing debug messages

#vcf_CXXFLAGS= -std=c++0x -g -pg -O3 -Wall -DNDEBUG   #options for not showing debug messages
#vcf_CXXFLAGS= -std=c++0x -g -O0 -Wall -DNDEBUG   #options for not showing debug messages


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

debug_src = scrm/random/constant_generator.cc scrm/random/constant_generator.h \
            scrm/forest-debug.cc \
            debug/particle-debug.cpp \
            debug/count-debug.cpp \
            debug/vcf-debug.cpp


arg_src = coalevent.cpp \
          particle.cpp \
          particleContainer.cpp \
          count.cpp count-semi-debug.cpp
       
param_src = pattern.cpp pfparam.cpp vcf.cpp 

#scrm_dummy_SOURCES = $(scrm_param) $(scrm_src) scrm/scrm.cc
#param_SOURCES = $(scrm_param) $(scrm_src)  $(param_src)  param_dummy.cpp
#coalevent_SOURCES = $(scrm_param) $(scrm_src)  coalevent_dummy.cpp
#particleCountainer_SOURCES = coalevent.cpp \
          #particle.cpp \
          #particleContainer.cpp $(scrm_param) $(scrm_src) $(param_src) \
          #paritcle_Container_dummy.cpp
#particle_SOURCES = coalevent.cpp \
          #particle.cpp $(scrm_param) $(scrm_src)  \
          #particle_dummy.cpp

shared_SOURCES = $(arg_src) $(scrm_param) $(scrm_src) $(param_src) 

pf_ARG_SOURCES = $(shared_SOURCES)  pf-ARG.cpp

pf_ARG_prof_CXXFLAGS = $(shared_SOURCES)  pf-ARG.cpp

pf_ARG_dbg_SOURCES = $(pf_ARG_SOURCES) $(debug_src)

unit_tests_SOURCES = tests/test_Vcf.cpp\
                     $(shared_SOURCES) \
                     tests/test_runner.cpp \
                     tests/test_param.cpp \
                     tests/test_pattern.cpp \
                     tests/test_random_generator.cpp \
                     tests/test_Coalevent.cpp \
                     tests/test_FroestState.cpp\
                     tests/test_PfCount.cpp\
                     tests/test_Vparticle.cpp\
                     #tests/test_Count.cpp\
                     #$(pf_ARG_SOURCES)

unit_tests_CXXFLAGS = -DUNITTEST -DNDEBUG -std=c++0x 

unit_tests_LDADD= -L/opt/local/lib -lcppunit -ldl  #link the cppunit unittest library in mac, cppunit was installed via macports
#unit_tests_LDADD= -L/usr/include 

include_HEADERS = vcf.hpp pfparam.hpp particle.hpp pattern.hpp\
                  scrm/forest.h scrm/model.h scrm/node.h scrm/node_container.h \
                  scrm/time_interval.h scrm/tree_point.h \
                  scrm/random/constant_generator.h scrm/random/mersenne_twister.h \
                  scrm/random/random_generator.h scrm/event.h \
                  scrm/summary_statistics/newick_tree.h \
                  scrm/summary_statistics/frequency_spectrum.h \
                  scrm/summary_statistics/seg_sites.h \
                  scrm/summary_statistics/tmrca.h
				  
