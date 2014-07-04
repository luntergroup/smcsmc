
#include <iostream>
#include <ctime>

#include "scrm/param.h"
#include "scrm/forest.h"
#include "scrm/random/random_generator.h"
#include "scrm/random/mersenne_twister.h"
#include <omp.h>


#ifndef UNITTEST
int main(int argc, char *argv[]){
  try {
    // Organize output
    std::ostream *output = &std::cout;

    Param user_para(argc, argv);

    Model top_model;
    user_para.parse(top_model);

    // Print help if user asked for it
    if (user_para.help()) {
      user_para.printHelp(*output); 
      return EXIT_SUCCESS;
    }
    if (user_para.version()) {
      *output << "scrm " << VERSION << std::endl;
      return EXIT_SUCCESS;
    }

    MersenneTwister rg = MersenneTwister(user_para.random_seed());

    //*output << "scrm " << VERSION << " |" << user_para << std::endl;
    *output << user_para << std::endl;
    *output << rg.seed() << std::endl;

    // Loop over the independent samples
    #pragma omp parallel for schedule(dynamic) 
    for (size_t rep_i=0; rep_i < top_model.loci_number(); ++rep_i) {

      Model* tmp_model = new Model(top_model);
      cout <<" doing " <<rep_i<<"th rep"<<endl;
      // Mark the start of a new independent sample
      //*output << std::endl << "//" << std::endl;

      // Now set up the ARG, and sample the initial tree
      Forest forest = Forest(tmp_model, &rg);
      forest.buildInitialTree();
      //forest.printSegmentSumStats(*output);

      while (forest.next_base() < top_model.loci_length()) { 
        // Sample next genealogy
        forest.sampleNextGenealogy();
        //forest.printSegmentSumStats(*output);
      }

      //forest.printLocusSumStats(*output);
    }

    // Clean-up and exit
    rg.clearFastFunc();
    return EXIT_SUCCESS;
  }
  catch (const exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    std::cerr << "Try 'scrm --help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
}
#endif
