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
    CPPUNIT_TEST(testWriteLog);
    CPPUNIT_TEST(testInsertMutationRateInScrmInput);
    CPPUNIT_TEST(testInsertSampleSizeInScrmInput);
    CPPUNIT_TEST(testInsertRecombRateAndSeqlenInScrmInput);
    CPPUNIT_TEST(testParse);
    CPPUNIT_TEST(testPrintVersion);
    CPPUNIT_TEST(testPrintHelp);
    CPPUNIT_TEST(testOutOfRange);
    CPPUNIT_TEST(testWrongType);
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
        CPPUNIT_ASSERT_EQUAL(this->input_->calibrate_lag, false);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(this->input_->lag_fraction, 0.0, this->epsilon_);
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
        CPPUNIT_ASSERT(this->input_->input_SegmentDataFileName == "");
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


    void testWriteLog(){
        CPPUNIT_ASSERT_EQUAL ( this->input_->log(), (int)0 );
    }


    void testInsertMutationRateInScrmInput(){
        // Default
        CPPUNIT_ASSERT_NO_THROW(this->input_->insertMutationRateInScrmInput());
        CPPUNIT_ASSERT(this->input_->scrm_input == "-t 8000.000000 "); // Note, there is a space at the end!

        // Parsed
        CPPUNIT_ASSERT_NO_THROW(this->input_->reInit());
        this->input_->scrm_input = "-t 800";
        CPPUNIT_ASSERT_NO_THROW(this->input_->insertMutationRateInScrmInput());
        CPPUNIT_ASSERT(this->input_->scrm_input == "-t 800");
    }


    void testInsertSampleSizeInScrmInput(){
        // Default
        CPPUNIT_ASSERT_NO_THROW(this->input_->insertSampleSizeInScrmInput());
        CPPUNIT_ASSERT(this->input_->scrm_input == "2 1 "); // Note, there is a space at the end!

        // Parsed
        CPPUNIT_ASSERT_NO_THROW(this->input_->reInit());
        this->input_->default_nsam = 10;
        CPPUNIT_ASSERT_NO_THROW(this->input_->insertSampleSizeInScrmInput());
        CPPUNIT_ASSERT(this->input_->scrm_input == "10 1 ");
    }


    void testInsertRecombRateAndSeqlenInScrmInput(){
        // Default
        CPPUNIT_ASSERT_NO_THROW(this->input_->insertRecombRateAndSeqlenInScrmInput());
        CPPUNIT_ASSERT(this->input_->scrm_input == "-r 800.000000 20000000 "); // Note, there is a space at the end!

        // Parsed
        CPPUNIT_ASSERT_NO_THROW(this->input_->reInit());
        this->input_->scrm_input = "-r 100 1000000";
        CPPUNIT_ASSERT_NO_THROW(this->input_->insertRecombRateAndSeqlenInScrmInput());
        CPPUNIT_ASSERT(this->input_->scrm_input == "-r 100 1000000");
        CPPUNIT_ASSERT_DOUBLES_EQUAL(this->input_->default_loci_length, 1000000, this->epsilon_);
    }


    void testParse(){
        char *argv1[] = { "./smcsmc",
                         "-Np", "1000",
                         "-nsam", "5",
                         "-ESS", "0.1",
                         "-EM", "5",
                         "-cap", "2000",
                         "-o", "tmp1",
                         "-tmax", "5",
                         "-p", "3*1+1",
                         "-lag", "0.5",
                         "-calibrate_lag", "0.3",
                         "-log",
                         "-heat",
                         "-seg", "noSuchFile",
                         "-unKnown"};
        // unKnown
        CPPUNIT_ASSERT_THROW(this->input_->parse(26, argv1), std::invalid_argument);

        // seg
        CPPUNIT_ASSERT_THROW(this->input_->parse(25, argv1), std::invalid_argument);
        CPPUNIT_ASSERT(this->input_->input_SegmentDataFileName == "noSuchFile");
        CPPUNIT_ASSERT_THROW(this->input_->parse(24, argv1), NotEnoughArg);

        // heat
        CPPUNIT_ASSERT_NO_THROW(this->input_->parse(23, argv1));
        CPPUNIT_ASSERT_EQUAL(this->input_->heat_bool, true);

        // log
        CPPUNIT_ASSERT_NO_THROW(this->input_->parse(22, argv1));
        CPPUNIT_ASSERT_EQUAL(this->input_->log_bool, true);

        // calibrate_lag
        CPPUNIT_ASSERT_NO_THROW(this->input_->parse(21, argv1));
        CPPUNIT_ASSERT_DOUBLES_EQUAL(this->input_->lag_fraction, 0.3, this->epsilon_);
        CPPUNIT_ASSERT_EQUAL(this->input_->calibrate_lag, true);
        CPPUNIT_ASSERT_THROW(this->input_->parse(20, argv1), NotEnoughArg);

        // lag
        CPPUNIT_ASSERT_NO_THROW(this->input_->parse(19, argv1));
        CPPUNIT_ASSERT_DOUBLES_EQUAL(this->input_->lag, 0.5, this->epsilon_);
        CPPUNIT_ASSERT_THROW(this->input_->parse(18, argv1), NotEnoughArg);

        // pattern
        CPPUNIT_ASSERT_NO_THROW(this->input_->parse(17, argv1));
        CPPUNIT_ASSERT(this->input_->pattern == "3*1+1");
        CPPUNIT_ASSERT_THROW(this->input_->parse(16, argv1), NotEnoughArg);

        // tmax
        CPPUNIT_ASSERT_NO_THROW(this->input_->parse(15, argv1));
        CPPUNIT_ASSERT_DOUBLES_EQUAL(this->input_->top_t(), 5, this->epsilon_);
        CPPUNIT_ASSERT_THROW(this->input_->parse(14, argv1), NotEnoughArg);

        // prefix
        CPPUNIT_ASSERT_NO_THROW(this->input_->parse(13, argv1));
        CPPUNIT_ASSERT(this->input_->out_NAME_prefix == "tmp1");
        CPPUNIT_ASSERT_THROW(this->input_->parse(12, argv1), NotEnoughArg);

        // cap
        CPPUNIT_ASSERT_NO_THROW(this->input_->parse(11, argv1));
        CPPUNIT_ASSERT_DOUBLES_EQUAL(this->input_->Ne_cap, 2000.0, this->epsilon_);
        CPPUNIT_ASSERT_EQUAL(this->input_->useCap, true);
        CPPUNIT_ASSERT_THROW(this->input_->parse(10, argv1), NotEnoughArg);

        // EM
        CPPUNIT_ASSERT_NO_THROW(this->input_->parse(9, argv1));
        CPPUNIT_ASSERT_EQUAL(this->input_->EM_steps, (int)5);
        CPPUNIT_ASSERT_THROW(this->input_->parse(8, argv1), NotEnoughArg);

        // ESS
        CPPUNIT_ASSERT_NO_THROW(this->input_->parse(7, argv1));
        CPPUNIT_ASSERT_EQUAL(this->input_->ESS(), 0.1);
        CPPUNIT_ASSERT_THROW(this->input_->parse(6, argv1), NotEnoughArg);

        // nsam
        CPPUNIT_ASSERT_NO_THROW(this->input_->parse(5, argv1));
        CPPUNIT_ASSERT_EQUAL(this->input_->default_nsam, (size_t)5);
        CPPUNIT_ASSERT_THROW(this->input_->parse(4, argv1), NotEnoughArg);

        // Np
        CPPUNIT_ASSERT_NO_THROW(this->input_->parse(3, argv1));
        CPPUNIT_ASSERT_EQUAL(this->input_->N, (size_t)1000);
        CPPUNIT_ASSERT_THROW(this->input_->parse(2, argv1), NotEnoughArg);
    }

    void testOutOfRange(){
        char *argv1[] = { "./smcsmc",
                         "-ESS", "2.0"};
        CPPUNIT_ASSERT_THROW(this->input_->parse(3, argv1), OutOfRange);

        char *argv3[] = { "./smcsmc",
                         "-calibrate_lag", "2.0"};
        CPPUNIT_ASSERT_THROW(this->input_->parse(3, argv3), OutOfRange);
    }

    void testWrongType(){
        char *argv1[] = { "./smcsmc",
                         "-Np", "asdf"};
        CPPUNIT_ASSERT_THROW(this->input_->parse(3, argv1), WrongType);
    }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION(TestPfParam);
