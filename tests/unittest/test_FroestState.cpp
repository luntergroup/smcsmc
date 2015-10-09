#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "particle.hpp"
#include "mersenne_twister.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestForestState : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestForestState );
  CPPUNIT_TEST ( testInitialization );
  CPPUNIT_TEST ( testCreateExampleTree );
  CPPUNIT_TEST ( testIncludeHaplotypesAtTips);
  CPPUNIT_TEST ( testTrackLocalTreeBranchLength );
  //CPPUNIT_TEST ( testLikelihood );
  CPPUNIT_TEST_SUITE_END();

 private:
    MersenneTwister *topRg;
    vector<int> record_event_in_epoch;

 public:
    void setUp() {
        topRg = new MersenneTwister(1234);
    }

    void tearDown() {
    }

    void testInitialization(){
        MersenneTwister *testRg = new MersenneTwister();
        ForestState testState = ForestState(new Model(4), testRg, record_event_in_epoch, false);
        CPPUNIT_ASSERT( testState.model().sample_size() == 4 );
        CPPUNIT_ASSERT( testState.random_generator() == testRg );
        CPPUNIT_ASSERT_EQUAL(0.0, testState.site_where_weight_was_updated() );
        CPPUNIT_ASSERT_EQUAL(1.0, testState.weight());
    }

    void testCreateExampleTree(){
        ForestState testState = ForestState(new Model(4), new MersenneTwister(1), record_event_in_epoch, true);
        CPPUNIT_ASSERT_NO_THROW( testState.createExampleTree());
        //cout<<endl;
        //CPPUNIT_ASSERT_NO_THROW(testState.printTree_cout());
        CPPUNIT_ASSERT_EQUAL(size_t(9), testState.nodes()->size());
        CPPUNIT_ASSERT( testState.local_root() == testState.nodes()->get(8) );
        CPPUNIT_ASSERT( testState.primary_root() == testState.nodes()->get(8) );
        CPPUNIT_ASSERT( testState.getLocalTreeLength() == 24 );
        CPPUNIT_ASSERT( testState.checkTree() == 1 );
    }

    void testIncludeHaplotypesAtTips(){
        ForestState testState = ForestState(new Model(4), new MersenneTwister(1), record_event_in_epoch, true);
        CPPUNIT_ASSERT_NO_THROW( testState.createExampleTree());
        int x[] = {1,0,0,1};
        vector <int> haplotypesAtTips(begin(x), end(x));
        CPPUNIT_ASSERT_NO_THROW( testState.include_haplotypes_at_tips(haplotypesAtTips) );
        for ( auto i = 0; i < haplotypesAtTips.size(); i++){
            CPPUNIT_ASSERT_EQUAL(haplotypesAtTips[i], testState.nodes()->get(i)->mutation_state());
        }
        // {1, -1, -1, 1}
        haplotypesAtTips[1] = -1;
        haplotypesAtTips[2] = -1;
        CPPUNIT_ASSERT_NO_THROW( testState.include_haplotypes_at_tips(haplotypesAtTips) );
        for ( auto i = 0; i < haplotypesAtTips.size(); i++){
            CPPUNIT_ASSERT_EQUAL(haplotypesAtTips[i], testState.nodes()->get(i)->mutation_state());
        }
    }

    void testTrackLocalTreeBranchLength(){
        // Need to check the local branch length is the following case:
        // (1). 0
        // (2). 2
        // (3). 6
        // (4). 20
        // (5). 21
        // (6). 24
        ForestState testState = ForestState(new Model(4), new MersenneTwister(1),record_event_in_epoch, true);
        CPPUNIT_ASSERT_NO_THROW( testState.createExampleTree());
        vector <int> haplotypesAtTips(4,1); // [1,1,1,1]
        CPPUNIT_ASSERT_NO_THROW( testState.include_haplotypes_at_tips(haplotypesAtTips) );
        CPPUNIT_ASSERT_EQUAL( (double)24 , testState.trackLocalTreeBranchLength());

        haplotypesAtTips[0] = 0; // [0,1,1,1]
        testState.include_haplotypes_at_tips(haplotypesAtTips);
        CPPUNIT_ASSERT_EQUAL( (double)24 , testState.trackLocalTreeBranchLength());

        haplotypesAtTips[0] = -1; // [-1,1,1,1]
        testState.include_haplotypes_at_tips(haplotypesAtTips);
        cout<<endl;
        CPPUNIT_ASSERT_NO_THROW(testState.printTree_cout());

        CPPUNIT_ASSERT_EQUAL( (double)21 , testState.trackLocalTreeBranchLength());

        haplotypesAtTips[0] = 0; // [0,1,1,1]
        haplotypesAtTips[2] = 0; // [0,1,0,1]
        testState.include_haplotypes_at_tips(haplotypesAtTips);
        CPPUNIT_ASSERT_EQUAL( (double)24 , testState.trackLocalTreeBranchLength());

        haplotypesAtTips[2] = -1; // [0,1,-1,1]
        testState.include_haplotypes_at_tips(haplotypesAtTips);
        CPPUNIT_ASSERT_EQUAL( (double)23 , testState.trackLocalTreeBranchLength());

        haplotypesAtTips[3] = -1; // [0,1,-1,-1]
        testState.include_haplotypes_at_tips(haplotypesAtTips);
        CPPUNIT_ASSERT_EQUAL( (double)6 , testState.trackLocalTreeBranchLength());

        haplotypesAtTips[0] = -1;
        haplotypesAtTips[1] = -1;
        haplotypesAtTips[2] = 0;
        haplotypesAtTips[3] = 1; // [-1,-1,0,1]
        testState.include_haplotypes_at_tips(haplotypesAtTips);
        CPPUNIT_ASSERT_EQUAL( (double)2, testState.trackLocalTreeBranchLength());

        haplotypesAtTips[3] = -1; // [-1,-1,0,-1]
        testState.include_haplotypes_at_tips(haplotypesAtTips);
        CPPUNIT_ASSERT_EQUAL( (double)0, testState.trackLocalTreeBranchLength());

        haplotypesAtTips[1] = 1; // [-1,1,0,-1]
        testState.include_haplotypes_at_tips(haplotypesAtTips);
        CPPUNIT_ASSERT_EQUAL( (double)20, testState.trackLocalTreeBranchLength());

        haplotypesAtTips[1] = -1;
        haplotypesAtTips[2] = -1; // [-1,-1,-1,-1]
        testState.include_haplotypes_at_tips(haplotypesAtTips);
        CPPUNIT_ASSERT_EQUAL( (double)0, testState.trackLocalTreeBranchLength());

        haplotypesAtTips[3] = 0;
        haplotypesAtTips[0] = 1; // [1,-1,-1,0]
        testState.include_haplotypes_at_tips(haplotypesAtTips);
        CPPUNIT_ASSERT_EQUAL( (double)20, testState.trackLocalTreeBranchLength());

    }


    void testLikelihood(){
        ForestState testState = ForestState(new Model(4), new MersenneTwister(1), record_event_in_epoch, true);
        CPPUNIT_ASSERT_NO_THROW( testState.createExampleTree());
        cout<<endl;
        vector <int> haplotypesAtTips(4,1); // [1,1,1,1]
        testState.include_haplotypes_at_tips(haplotypesAtTips);
        CPPUNIT_ASSERT_NO_THROW(testState.printTree_cout());
        haplotypesAtTips[0] = 0; // [0,1,1,1]
        testState.include_haplotypes_at_tips(haplotypesAtTips);
        CPPUNIT_ASSERT_NO_THROW(testState.printTree_cout());
        haplotypesAtTips[2] = 0; // [0,1,0,1]
        testState.include_haplotypes_at_tips(haplotypesAtTips);
        CPPUNIT_ASSERT_NO_THROW(testState.printTree_cout());
    }


};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestForestState );

