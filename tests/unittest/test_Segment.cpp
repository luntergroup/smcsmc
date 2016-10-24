#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "segdata.hpp"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestSegment : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE(TestSegment);
    CPPUNIT_TEST(testMainConstructor);
    CPPUNIT_TEST(testWrongNumberOfEntry);
    CPPUNIT_TEST(testInvalidSegmentStartPosition);
    CPPUNIT_TEST(testInvalidInputFile);
    CPPUNIT_TEST_SUITE_END();

  private:
    Segment* segFile_;
    double epsilon_;

  public:
    void setUp() {
        this->segFile_ = new Segment ();
        this->epsilon_ = 1e-10;
    }

    void tearDown() {
        delete segFile_;
    }


    void testMainConstructor(){
        CPPUNIT_ASSERT_NO_THROW(Segment( "data/exampleData/example.seg", (size_t)4, 1000000.0, 100.0));
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->init( "data/exampleData/example.seg", (size_t)4, 1000000.0, 100.0));
        //1       521     T       F       1       01.0
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->read_new_line());
        //522     2721    T       F       1       0111
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->read_new_line());
        //3243    1758    T       F       1       10//
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->read_new_line());
        //5001    1296    T       F       1       0000
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->read_new_line());
        //6297    1       T       F       1       ....
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->read_new_line());
        //6298    4669    T       F       1       0110
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->read_new_line());
        //10967   880     T       T       1       0100
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->read_new_line());
        //1       708     T       F       2       1010
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->read_new_line());
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->read_new_line());

        CPPUNIT_ASSERT_NO_THROW(this->segFile_->init( "", (size_t)4, 1000000.0, 100.0));
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->read_new_line());
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->read_new_line());
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->read_new_line());
    }


    void testWrongNumberOfEntry(){
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->init( "data/exampleData/example.WrongNumberOfEntry.seg", (size_t)4, 1000000.0, 100.0));
        //1       521     T       F       1       01.0
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->read_new_line());
        //522     2721    T       F       1       01111
        CPPUNIT_ASSERT_THROW(this->segFile_->read_new_line(), WrongNumberOfEntry);
    }


    void testInvalidSegmentStartPosition(){
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->init( "data/exampleData/example.InvalidSegmentStartPosition.seg", (size_t)4, 1000000.0, 100.0));
        //1	521	T	F	1	01.0
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->read_new_line());
        //522	2721	T	F	1	0111
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->read_new_line());
        //3243	1758	T	F	1	10
        CPPUNIT_ASSERT_NO_THROW(this->segFile_->read_new_line());
        //5003	1296	T	F	1	0000
        CPPUNIT_ASSERT_THROW(this->segFile_->read_new_line(), InvalidSegmentStartPosition);
    }


    void testInvalidInputFile(){
        CPPUNIT_ASSERT_THROW(this->segFile_->init("noSuchFile", (size_t)4, 1000000.0, 100.0), InvalidInputFile);
        CPPUNIT_ASSERT_THROW(Segment("noSuchFile", (size_t)4, 1000000.0, 100.0), InvalidInputFile);
    }

};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION(TestSegment);
