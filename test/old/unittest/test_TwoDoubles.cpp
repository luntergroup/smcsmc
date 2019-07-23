#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "count.hpp"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestTwoDoubles : public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE( TestTwoDoubles );
    CPPUNIT_TEST( testMainConstructor );
    CPPUNIT_TEST( testAdd );
    CPPUNIT_TEST_SUITE_END();

 private:
    TwoDoubles * testNumber1_;
    TwoDoubles * testNumber2_;
    double epsilon_8, epsilon_7;

 public:
    void setUp() {
        this->testNumber1_ = new TwoDoubles();
        this->testNumber2_ = new TwoDoubles(3.1);
        this->epsilon_8 = 1e-8;
        this->epsilon_7 = 1e-7;
    }

    void tearDown() {
        delete testNumber1_;
        delete testNumber2_;
    }

    void testMainConstructor(){
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.0, this->testNumber1_->small_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.0, this->testNumber1_->big_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 3.1, this->testNumber2_->small_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.0, this->testNumber2_->big_, this->epsilon_8 );
    }

    void testAdd() {
        CPPUNIT_ASSERT_NO_THROW (this->testNumber1_->add(0.4) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.0, this->testNumber1_->small_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.4, this->testNumber1_->big_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.4, this->testNumber1_->final_answer(), this->epsilon_8 );

        CPPUNIT_ASSERT_NO_THROW (this->testNumber1_->add(67108864) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.0, this->testNumber1_->small_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 67108864.4, this->testNumber1_->big_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 67108864.4, this->testNumber1_->final_answer(), this->epsilon_8 );

        CPPUNIT_ASSERT_NO_THROW (this->testNumber1_->add(0.1) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.1, this->testNumber1_->small_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 67108864.4, this->testNumber1_->big_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 67108864.5, this->testNumber1_->final_answer(), this->epsilon_8 );

        CPPUNIT_ASSERT_NO_THROW (this->testNumber1_->add(0.001) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.101, this->testNumber1_->small_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 67108864.4, this->testNumber1_->big_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 67108864.501, this->testNumber1_->final_answer(), this->epsilon_8 );

        CPPUNIT_ASSERT_NO_THROW (this->testNumber1_->add(0.898) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.999, this->testNumber1_->small_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 67108864.4, this->testNumber1_->big_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 67108865.399, this->testNumber1_->final_answer(), this->epsilon_8 );

        // Starting to lose accuracy
        CPPUNIT_ASSERT_NO_THROW (this->testNumber1_->add(0.002) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.0, this->testNumber1_->small_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 67108865.401, this->testNumber1_->big_, this->epsilon_7 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 67108865.401, this->testNumber1_->final_answer(), this->epsilon_7 );

        CPPUNIT_ASSERT_NO_THROW (this->testNumber1_->add(0.002) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.002, this->testNumber1_->small_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 67108865.401, this->testNumber1_->big_, this->epsilon_7 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 67108865.403, this->testNumber1_->final_answer(), this->epsilon_7 );

        CPPUNIT_ASSERT_NO_THROW (this->testNumber2_->add(0.4) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.0, this->testNumber2_->small_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 3.5, this->testNumber2_->big_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 3.5, this->testNumber2_->final_answer(), this->epsilon_8 );

        CPPUNIT_ASSERT_NO_THROW (this->testNumber2_->add(67108864) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.0, this->testNumber2_->small_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 67108867.5, this->testNumber2_->big_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 67108867.5, this->testNumber2_->final_answer(), this->epsilon_8 );

        // Starting to lose accuracy
        CPPUNIT_ASSERT_NO_THROW (this->testNumber2_->add(0.001) );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0.001, this->testNumber2_->small_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 67108867.5, this->testNumber2_->big_, this->epsilon_8 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL ( 67108867.501, this->testNumber2_->final_answer(), this->epsilon_7 );
    }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestTwoDoubles );
