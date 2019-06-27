#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "particle.hpp"
#include "mersenne_twister.h"
#include "constant_generator.h"
#include "forest.h"
//#include "../../src/random/constant_generator.h"
//#include "../../src/random/mersenne_twister.h"
#include "event.h"
//#include "../../src/summary_statistics/tmrca.h"


#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestForestState : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestForestState );
  CPPUNIT_TEST ( testInitialization );
  CPPUNIT_TEST ( testCreateExampleTree );
  //CPPUNIT_TEST ( testIncludeHaplotypesAtTips);
  //CPPUNIT_TEST ( testTrackLocalTreeBranchLength );
  //CPPUNIT_TEST ( testCopyExampleTree );
  //CPPUNIT_TEST ( testLikelihood );
  CPPUNIT_TEST_SUITE_END();

 private:
    RandomGenerator *topRg;
    vector<int> record_event_in_epoch;
  Model *model, *model_2pop;
  Forest *forest, *forest_2pop;
  MersenneTwister *rg;

 public:
    void setUp() {
        //topRg = new RandomGenerator();
        topRg = new MersenneTwister(true, 1234);
            //rg = new MersenneTwister(1234);
    //model = new Model(5);
    //model_2pop = new Model(5);
    //model_2pop->set_population_number(2);
    //forest = new Forest(model, rg);
    //forest->createExampleTree();
    //forest_2pop = new Forest(model_2pop, rg);
    //forest_2pop->createExampleTree();

    }

    void tearDown() {
        delete topRg;
            //delete forest, forest_2pop;
    //delete model, model_2pop;
    //delete rg;

    }

    void testInitialization(){
        //RandomGenerator *testRg = new MersenneTwister(true, 1);
        Model *testModel = new Model(4);
        ForestState testState = ForestState(testModel, topRg, record_event_in_epoch, false);
        CPPUNIT_ASSERT( testState.model().sample_size() == 4 );
        CPPUNIT_ASSERT( testState.random_generator() == topRg );
        CPPUNIT_ASSERT_EQUAL(0.0, testState.site_where_weight_was_updated() );
        CPPUNIT_ASSERT_EQUAL(1.0, testState.weight());
        delete testModel;
        //delete testRg;
    }

    void testCreateExampleTree(){
        Model *testModel = new Model(4);
        //MersenneTwister *testRg = new MersenneTwister(1);
        ForestState testState = ForestState(testModel, topRg, record_event_in_epoch, false);
        CPPUNIT_ASSERT_NO_THROW( testState.createExampleTree());
        //cout<<endl;
        //CPPUNIT_ASSERT_NO_THROW(testState.printTree_cout());
        CPPUNIT_ASSERT_EQUAL(size_t(9), testState.nodes()->size());
        CPPUNIT_ASSERT( testState.local_root() == testState.nodes()->get(8) );
        CPPUNIT_ASSERT( testState.primary_root() == testState.nodes()->get(8) );
        CPPUNIT_ASSERT( testState.getLocalTreeLength() == 24 );
        CPPUNIT_ASSERT( testState.checkTree() == 1 );
        delete testModel;
        //delete testRg;
    }

    void testCopyExampleTree(){
        Model *testModel = new Model(4);
        MersenneTwister *testRg = new MersenneTwister(1);
        ForestState testState = ForestState(testModel, testRg, record_event_in_epoch, true);
        CPPUNIT_ASSERT_NO_THROW( testState.createExampleTree());
        //ForestState * copyedTree = new ForestState(testState);
        //cout<<endl;
        //CPPUNIT_ASSERT_NO_THROW(testState.printTree_cout());
        //CPPUNIT_ASSERT_EQUAL(size_t(9), testState.nodes()->size());
        //CPPUNIT_ASSERT( testState.local_root() == testState.nodes()->get(8) );
        //CPPUNIT_ASSERT( testState.primary_root() == testState.nodes()->get(8) );
        //CPPUNIT_ASSERT( testState.getLocalTreeLength() == 24 );
        //CPPUNIT_ASSERT( testState.checkTree() == 1 );
        delete testModel;
        delete testRg;
    }

    void testIncludeHaplotypesAtTips(){
        Model *testModel = new Model(4);
        MersenneTwister *testRg = new MersenneTwister(1);
        ForestState testState = ForestState(testModel, testRg, record_event_in_epoch, true);
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
        delete testModel;
        delete testRg;
    }

    void testTrackLocalTreeBranchLength(){
        // Need to check the local branch length is the following case:
        // (1). 0
        // (2). 2
        // (3). 6
        // (4). 20
        // (5). 21
        // (6). 24
        Model *testModel = new Model(4);
        MersenneTwister *testRg = new MersenneTwister(1);
        ForestState testState = ForestState(testModel, testRg, record_event_in_epoch, true);
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

        delete testModel;
        delete testRg;

    }


    void testLikelihood(){
        Model *testModel = new Model(4);
        MersenneTwister *testRg = new MersenneTwister(1);
        ForestState testState = ForestState(testModel, testRg, record_event_in_epoch, false);
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
        delete testModel;
        delete testRg;
    }


};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestForestState );

