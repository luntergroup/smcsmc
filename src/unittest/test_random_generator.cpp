/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "random/mersenne_twister.h"
#include "random/random_generator.h"

class TestRandomGenerator : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestRandomGenerator );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testSampleExpoLimit );
  CPPUNIT_TEST( testSampleUnitExpo );
  CPPUNIT_TEST( testSampleExpo );
  CPPUNIT_TEST( testSampleExpoExpoLimit );
  CPPUNIT_TEST( testSampleInt );
  CPPUNIT_TEST( testSeeding );

  CPPUNIT_TEST_SUITE_END();

 private:
  RandomGenerator *rg;

 public:
  void setUp() {
    rg = new MersenneTwister(5);
  }

  void tearDown() {
    delete rg;
  }

  void testConstructor() {
    CPPUNIT_ASSERT_EQUAL( (size_t)5, rg->seed() );
    CPPUNIT_ASSERT( rg->unit_exponential_ > 0 );
  }

  void testSampleExpoLimit() {
    double expo;
    for (size_t i = 0; i < 1000; ++i) {
      expo = rg->sampleExpoLimit(1, 2);
      CPPUNIT_ASSERT( expo == -1 || expo > 0 );
      CPPUNIT_ASSERT( expo <= 2 ); 
    }
  }

  void testSampleUnitExpo() {
    size_t n = 10000;
    double expo = 0.0;
    for (size_t i = 0; i < n; ++i) {
      expo += rg->sampleUnitExponential();
    }
    expo /= n;
    CPPUNIT_ASSERT( 0.99 <= expo && expo <= 1.01 );
  }

  void testSampleExpo() {
    size_t n = 10000;
    double expo = 0;
    for (size_t i = 0; i < n; ++i) {
      expo += rg->sampleExpo(5);
    }
    expo /= n;
    CPPUNIT_ASSERT( 0.199 <= expo && expo <= 0.201 );
  }
  
  void testSampleExpoExpoLimit() {
    size_t n = 10000;
    double expo = 0, sample = 0;
    size_t sample_number = 0;
    
    // c = 0;
    for (size_t i = 0; i < n; ++i) {
      sample = rg->sampleExpoExpoLimit(2,0,1);
      if (sample >= 0) {
        expo += sample;
        ++sample_number;
      }
    }
    expo /= sample_number;
    // Expected: 0.34 
    CPPUNIT_ASSERT( 0.32 < expo && expo < 0.36 );

    expo = 0;
    sample_number = 0;
    for (size_t i = 0; i < n; ++i) {
      sample = rg->sampleExpoExpoLimit(2,0,2);
      if (sample >= 0) {
        expo += sample;
        ++sample_number;
      }
    }
    expo /= sample_number;
    // Expected: 0.46 
    CPPUNIT_ASSERT( 0.44 < expo && expo < 0.48 );

    expo = 0;
    sample_number = 0;
    for (size_t i = 0; i < n; ++i) {
      sample = rg->sampleExpoExpoLimit(2,0,4);
      if (sample >= 0) {
        expo += sample;
        ++sample_number;
      }
    }
    expo /= sample_number;
    // Expected: 0.50 
    CPPUNIT_ASSERT( 0.48 < expo && expo < 0.52 );
  }

  void testSampleInt() {
    size_t n = 50000;
    int sample;
    int result[5] = { 0 };

    for (size_t i = 0; i < n; ++i) {
      sample = rg->sampleInt(5);
      ++result[sample]; 
    }

    for (size_t i = 0; i < 5; ++i) {
      //std::cout << i << " : " << result[i] << std::endl;
      CPPUNIT_ASSERT( 9900 < result[i] && result[i] < 10100 ); 
    }
  }

  void testSeeding() {
    rg->set_seed(5);
    int sample = rg->sampleInt(10000);
    rg->set_seed(5);
    CPPUNIT_ASSERT_EQUAL( sample, rg->sampleInt(10000) );

    MersenneTwister rg2 = MersenneTwister(5);
    CPPUNIT_ASSERT_EQUAL( sample, rg2.sampleInt(10000) );
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestRandomGenerator );
