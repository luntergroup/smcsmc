#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "pattern.hpp"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestPattern : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE( TestPattern );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST( testExtractNumberOfSegment );
    CPPUNIT_TEST( testFaultyPattern );
    CPPUNIT_TEST_SUITE_END();

  private:
    Pattern * pattern_tester;
    string pattern;
    double top_t ;
    vector <double> t_i;
    vector <double> new_ti;

    void testExtractNumberOfSegmentCore(string pattern, size_t expectedNumSeg, size_t expectedNumGroup ){
        CPPUNIT_ASSERT_NO_THROW ( this->pattern_tester->extractFromPattern(pattern));
        CPPUNIT_ASSERT_EQUAL( expectedNumSeg, pattern_tester->num_seg_ );
        if ( expectedNumSeg > 1 ){
          t_i = pattern_tester->extract_Segment( );
          CPPUNIT_ASSERT_EQUAL(pattern_tester->top_t_, t_i.back());
          CPPUNIT_ASSERT_EQUAL(t_i.size(), pattern_tester->num_seg_);
          new_ti = pattern_tester->regroup_Segment ( t_i );
          CPPUNIT_ASSERT_EQUAL(expectedNumGroup, new_ti.size());
        }
        CPPUNIT_ASSERT_NO_THROW ( this->pattern_tester->init());
    }

  public:
    void setUp() {
        this->pattern_tester = new Pattern();
        top_t = 2.0;
    }


    void tearDown() {
        delete pattern_tester;
    }


    void testMainConstructor(){
        CPPUNIT_ASSERT_EQUAL(this->pattern_tester->pattern_str.size(), (size_t)0);
        CPPUNIT_ASSERT_EQUAL(this->pattern_tester->seg_level1_vec_.size(), (size_t)0);
        CPPUNIT_ASSERT_EQUAL(this->pattern_tester->seg_level2_vec_.size(), (size_t)0);
        CPPUNIT_ASSERT_EQUAL(this->pattern_tester->top_t_, 0.0);
        CPPUNIT_ASSERT(this->pattern_tester->expr_ == NULL);
        CPPUNIT_ASSERT_EQUAL(this->pattern_tester->num_seg_, (size_t)0);
        CPPUNIT_ASSERT_NO_THROW(Pattern("3+2+2+2+2+2+3", top_t));
    }


    void testFaultyPattern () {
        CPPUNIT_ASSERT_THROW( Pattern("1*", top_t), PatternDigitsExpected ); // Expect digits followed by *
        CPPUNIT_ASSERT_THROW( Pattern("1adf", top_t), PatternTimesExpected );
        CPPUNIT_ASSERT_THROW( Pattern("1adf*", top_t), PatternTimesExpected );
        CPPUNIT_ASSERT_THROW( Pattern("1*2a", top_t), PatternAddsExpected );
        CPPUNIT_ASSERT_THROW( Pattern("10-", top_t), PatternTimesExpected ); // Expect times followed by 10
        CPPUNIT_ASSERT_THROW( Pattern("-", top_t), PatternDigitsExpected );
        CPPUNIT_ASSERT_THROW( Pattern("+", top_t), PatternDigitsExpected );
        CPPUNIT_ASSERT_THROW( Pattern("-5", top_t), PatternDigitsExpected );
        CPPUNIT_ASSERT_THROW( Pattern("+4", top_t), PatternDigitsExpected );
    }


    void testExtractNumberOfSegment( ){
        this->testExtractNumberOfSegmentCore("3+2+2+2+2+2+3", (size_t)16, (size_t)7);
        this->testExtractNumberOfSegmentCore("1*4+25*2+1*4+1*6", (size_t)64, (size_t)28);
        this->testExtractNumberOfSegmentCore("0+0", (size_t)0, (size_t)2); // (size_t)2 is not checked ...
        this->testExtractNumberOfSegmentCore("0", (size_t)0, (size_t)1);
        this->testExtractNumberOfSegmentCore("1", (size_t)1, (size_t)1);
        this->testExtractNumberOfSegmentCore("3", (size_t)3, (size_t)1);
        this->testExtractNumberOfSegmentCore("1+1", (size_t)2, (size_t)2);
    }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestPattern );
