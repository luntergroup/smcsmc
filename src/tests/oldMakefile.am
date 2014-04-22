TESTS = unit_tests 

check_PROGRAMS = unit_tests 

src_dir = ../src/ 
# this sub does not work ...

scrm_dir = ../src/scrm/

scrm_param = ../src/scrm/param.cc ../src/scrm/param.h \
        ../src/scrm/model.cc ../src/scrm/model.h \
        ../src/scrm/summary_statistics/newick_tree.cc \
        ../src/scrm/summary_statistics/frequency_spectrum.cc \
        ../src/scrm/summary_statistics/seg_sites.cc \
		../src/scrm/summary_statistics/tmrca.cc \
		../src/scrm/summary_statistics/first_last_tmrca.cc


scrm_src = ../src/scrm/node.cc ../src/scrm/node.h \
        ../src/scrm/node_container.cc ../src/scrm/node_container.h \
        ../src/scrm/time_interval.cc ../src/scrm/time_interval.h \
        ../src/scrm/tree_point.cc ../src/scrm/tree_point.h \
        ../src/scrm/forest.cc \
        ../src/scrm/random/mersenne_twister.h  \
        ../src/scrm/random/random_generator.cc \
        ../src/scrm/random/mersenne_twister.cc   ../src/scrm/random/random_generator.h \
        ../src/scrm/random/fastfunc.cc           ../src/scrm/random/fastfunc.h 

debug_src = ../src/scrm/random/constant_generator.cc ../src/scrm/random/constant_generator.h \
            ../src/scrm/forest-debug.cc \
            ../src/particle-debug.cpp \
            ../src/count-debug.cpp \
            ../src/vcf-debug.cpp


arg_src = ../src/coalevent.hpp ../src/coalevent.cpp \
        ../src/particle.hpp ../src/particle.cpp \
        ../src/count.cpp ../src/count.hpp ../src/count-semi-debug.cpp 

vcf_src = ../src/vcf.hpp ../src/vcf.cpp 


param_src = ../src/pattern.cpp ../src/pfparam.cpp 

pf_ARG_SOURCES = $(debug_src) $(arg_src) $(scrm_param) $(scrm_src) $(vcf_src) $(param_src) 

unit_tests_SOURCES = test_Vcf.cpp\
                     $(pf_ARG_SOURCES) \
                     test_runner.cpp \
                     test_param.cpp \
                     test_random_generator.cpp
                     #test_Coalevent.cpp\
                     #test_FroestState.cpp\
                     #test_PfCount.cpp\
                     #test_Vparticle.cpp\
                     #test_Count.cpp\
                     #test_Vcf.cpp\
                     #$(pf_ARG_SOURCES)

unit_tests_CXXFLAGS = -DUNITTEST -DNDEBUG -std=c++0x 
#unit_tests_CXXFLAGS = -DUNITTEST  -std=c++0x 


unit_tests_LDADD= -L/opt/local/lib -lcppunit -ldl  #link the cppunit unittest library in mac, cppunit was installed via macports
#unit_tests_LDADD= -L/usr/include 

include_HEADERS = ../src/vcf.hpp ../src/pfparam.hpp ../src/particle.hpp \
                  ../src/scrm/forest.h ../src/scrm/model.h ../src/scrm/node.h ../src/scrm/node_container.h \
                  ../src/scrm/time_interval.h ../src/scrm/tree_point.h \
                  ../src/scrm/random/constant_generator.h ../src/scrm/random/mersenne_twister.h \
                  ../src/scrm/random/random_generator.h 
