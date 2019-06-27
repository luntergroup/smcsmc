#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

//#include "pfparam.hpp"
//#include "../src/usage.hpp"
//#include "variantReader.hpp"
//#include "param.h"
#include "model.h"
//#include "forest.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestPfCount : public CppUnit::TestCase {
  
  CPPUNIT_TEST_SUITE( TestPfCount );

  CPPUNIT_TEST( somemethod ); 
  
  CPPUNIT_TEST_SUITE_END();

 private:
    Model model;

 public:
  void somemethod() {
      } 
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestPfCount );
