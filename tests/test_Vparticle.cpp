#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <boost/lexical_cast.hpp> 

#include "../src/pfparam.hpp"
#include "../src/usage.hpp"
#include "../src/vcf.hpp"

#include "../src/scrm/param.h"
#include "../src/scrm/model.h"
#include "../src/scrm/forest.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestVparticle : public CppUnit::TestCase {
  
  CPPUNIT_TEST_SUITE( TestVparticle );

  CPPUNIT_TEST( somemethod ); 
  
  CPPUNIT_TEST_SUITE_END();

 private:
    Model model;

 public:
  void somemethod() {}
      
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestVparticle );
