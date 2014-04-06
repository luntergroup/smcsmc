bin_PROGRAMS = pf-ARG pf-ARG_dbg 
PROG = PF-ARG 


pf_ARG_CXXFLAGS= -std=c++0x -g -pg -O3 -Wall -DNDEBUG   #options for not showing debug messages
#pf_ARG_CXXFLAGS= -std=c++0x -g -O0 -Wall 

#pf_ARG_CXXFLAGS= -std=c++0x -g -O0 -Wall -DNDEBUG   #options for not showing debug messages
pf_ARG_dbg_CXXFLAGS= -std=c++0x  #options for showing debug messages

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
            particle-debug.cpp \
            count-debug.cpp \
            vcf-debug.cpp


arg_src = coalevent.hpp coalevent.cpp \
        particle.hpp particle.cpp \
        count.cpp count.hpp count-semi-debug.cpp\
        pf-ARG.cpp

vcf_src = vcf.hpp vcf.cpp 

param_src = pattern.cpp pfparam.cpp 

#usage_src = usage.hpp usage.cpp 

pf_ARG_SOURCES =  $(arg_src) $(scrm_param) $(scrm_src) $(vcf_src) $(param_src) $(debug_src)

pf_ARG_dbg_SOURCES = $(pf_ARG_SOURCES) 

#vcf_SOURCES = $(vcf_src) testvcf.cpp $(param_src) $(usage_src) $(scrm_param)

include_HEADERS = vcf.hpp pfparam.hpp particle.hpp pattern.hpp\
                  scrm/forest.h scrm/model.h scrm/node.h scrm/node_container.h \
                  scrm/time_interval.h scrm/tree_point.h \
                  scrm/random/constant_generator.h scrm/random/mersenne_twister.h \
                  scrm/random/random_generator.h scrm/event.h \
                  scrm/summary_statistics/newick_tree.h \
		  scrm/summary_statistics/frequency_spectrum.h \
		  scrm/summary_statistics/seg_sites.h \
		  scrm/summary_statistics/tmrca.h
				  
