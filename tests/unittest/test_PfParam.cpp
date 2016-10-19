#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "pfparam.hpp"
#include "segdata.hpp"
#include "param.h"
#include "model.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestPfParam : public CppUnit::TestCase {

    CPPUNIT_TEST_SUITE(TestPfParam);
    CPPUNIT_TEST(testMainConstructor);
    CPPUNIT_TEST(testReInit);
    CPPUNIT_TEST(testParse);
    CPPUNIT_TEST(testPrintVersion);
    CPPUNIT_TEST(testPrintHelp);
    CPPUNIT_TEST_SUITE_END();

 private:
    //Model model;
    PfParam* input_;
    double epsilon_;

  public:
    void setUp() {
        this->input_ = new PfParam ();
        this->epsilon_ = 1e-10;
    }

    void tearDown() {
        delete input_;
    }

    void testMainConstructor(){
        CPPUNIT_ASSERT_EQUAL(this->input_->default_nsam, (size_t)2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(this->input_->default_mut_rate, 1e-8, this->epsilon_);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(this->input_->default_recomb_rate, 1e-9, this->epsilon_);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(this->input_->default_loci_length, 2e7, this->epsilon_);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(this->input_->default_num_mut, this->input_->default_mut_rate*40000*this->input_->default_loci_length, this->epsilon_);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(this->input_->original_recombination_rate_, 0.0, this->epsilon_);
        CPPUNIT_ASSERT_EQUAL(this->input_->N, (size_t)100);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(this->input_->lag, 0.0, this->epsilon_);
        CPPUNIT_ASSERT(this->input_->out_NAME_prefix == "smcsmc");
        CPPUNIT_ASSERT_DOUBLES_EQUAL(this->input_->ESS(), 0.5, this->epsilon_);
        CPPUNIT_ASSERT_EQUAL(this->input_->ESS_default_bool, true);
        CPPUNIT_ASSERT_EQUAL(this->input_->log_bool, true);
        CPPUNIT_ASSERT_EQUAL(this->input_->heat_bool, false);
        CPPUNIT_ASSERT_EQUAL(this->input_->online_bool, false);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(this->input_->heat_seq_window, 400.0, this->epsilon_);
        CPPUNIT_ASSERT_EQUAL(this->input_->EM_steps, int(0));
        CPPUNIT_ASSERT_EQUAL(this->input_->EM_bool, false);
        CPPUNIT_ASSERT(this->input_->Segfile == NULL);
        CPPUNIT_ASSERT(this->input_->SCRMparam == NULL);
        CPPUNIT_ASSERT(this->input_->rg == NULL);
        CPPUNIT_ASSERT(this->input_->scrm_input == "");
        CPPUNIT_ASSERT_DOUBLES_EQUAL(this->input_->top_t_, 2.0, this->epsilon_);
        CPPUNIT_ASSERT_EQUAL(this->input_->EMcounter(), (size_t)0);
        CPPUNIT_ASSERT_EQUAL(this->input_->argc_, (int)0);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(this->input_->Ne_cap, 200000.0, this->epsilon_);
        CPPUNIT_ASSERT_EQUAL(this->input_->useCap, false);
        CPPUNIT_ASSERT_EQUAL(this->input_->help(), false);
        CPPUNIT_ASSERT_EQUAL(this->input_->version(), false);
    }

    void testReInit(){
        this->input_->Segfile = new Segment();
        this->input_->SCRMparam = new Param();
        this->input_->rg = new MersenneTwister();
        CPPUNIT_ASSERT_NO_THROW(this->input_->reInit());
        this->testMainConstructor();
    }



    void testPrintHelp(){
        char *argv[] = { "./smcsmc" };
        CPPUNIT_ASSERT_NO_THROW ( this->input_->parse(1, argv) );

        CPPUNIT_ASSERT_EQUAL((size_t)0, input_->argv_.size());
        CPPUNIT_ASSERT_EQUAL( true, input_->help());
        CPPUNIT_ASSERT_NO_THROW ( this->input_->printHelp() );

        char *argv1[] = { "./smcsmc", "-h" };
        CPPUNIT_ASSERT_NO_THROW ( this->input_->parse(2, argv1) );

        CPPUNIT_ASSERT_EQUAL((size_t)1, input_->argv_.size());
        CPPUNIT_ASSERT_EQUAL( true, input_->help());
        CPPUNIT_ASSERT_NO_THROW ( this->input_->printHelp() );

        char *argv2[] = { "./smcsmc", "-help" };
        CPPUNIT_ASSERT_NO_THROW ( this->input_->parse(2, argv2) );

        CPPUNIT_ASSERT_EQUAL((size_t)1, input_->argv_.size());
        CPPUNIT_ASSERT_EQUAL( true, input_->help());
        CPPUNIT_ASSERT_NO_THROW ( this->input_->printHelp() );
    }


    void testPrintVersion(){
        char *argv1[] = { "./smcsmc", "-v" };
        CPPUNIT_ASSERT_NO_THROW ( this->input_->parse(2, argv1) );
        CPPUNIT_ASSERT_EQUAL( true, input_->version());
        CPPUNIT_ASSERT_NO_THROW ( this->input_->printVersion(&std::cout) );

        char *argv2[] = { "./smcsmc", "-version" };
        CPPUNIT_ASSERT_NO_THROW ( this->input_->parse(2, argv2) );
        CPPUNIT_ASSERT_EQUAL( true, input_->version());
        CPPUNIT_ASSERT_NO_THROW ( this->input_->printVersion(&std::cout) );
    }


    void testParse(){

        char *argv1[] = { "./smcsmc",
                         "-ref", "data/testData/PG0390-C.test.ref",
                         "-alt", "data/testData/PG0390-C.test.alt",
                         "-plaf", "data/testData/labStrains.test.PLAF.txt",
                         "-panel", "data/testData/labStrains.test.panel.txt",
                         "-o", "tmp1",
                         "-exclude", "data/testData/labStrains.test.exclude.txt", "-vcfOut", "-z"};

        CPPUNIT_ASSERT_NO_THROW(this->input_->reInit());
    }

};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION(TestPfParam);
