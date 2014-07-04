#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <boost/lexical_cast.hpp> 

#include "../src/particle.hpp"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestForestState : public CppUnit::TestCase {
  
  CPPUNIT_TEST_SUITE( TestForestState );

  CPPUNIT_TEST( testInitialization );

  CPPUNIT_TEST( testCreateExampleTree );
  
    
  CPPUNIT_TEST_SUITE_END();

 private:
    //Model model;
    Model * model;
    MersenneTwister *rg;

 public:
 
    void testInitialization(){
        rg = new MersenneTwister();
        
        ForestState new_state = ForestState(new Model(4),rg);
        
        //CPPUNIT_ASSERT( new_state.model().sample_size() == 4 );
        //CPPUNIT_ASSERT( new_state.random_generator() == rg );
        //CPPUNIT_ASSERT_EQUAL(0.0, new_state.site_where_weight_was_updated() );
        //CPPUNIT_ASSERT_EQUAL(1.0, new_state.weight());
        //CPPUNIT_ASSERT(new_state.previous_state == NULL);
        cout<<endl;
        //CPPUNIT_ASSERT(new_state.printTree());
        
        //delete rg;
        //delete new_state.writable_model();
    }
     
    
    void testCreateExampleTree(){

        //ForestState new_state;
        //CPPUNIT_ASSERT_NO_THROW( new_state.createExampleTree());
        //cout<<endl;

        //CPPUNIT_ASSERT(new_state.printTree_cout());
        
        //CPPUNIT_ASSERT_EQUAL(size_t(9), new_state.nodes()->size());
        //CPPUNIT_ASSERT( new_state.local_root() == new_state.nodes()->get(8) );
        //CPPUNIT_ASSERT( new_state.primary_root() == new_state.nodes()->get(8) );
        //CPPUNIT_ASSERT( new_state.local_tree_length() == 24 );
        //CPPUNIT_ASSERT( new_state.checkTree() == 1 );
        //CPPUNIT_ASSERT_EQUAL(0.0, new_state.site_where_weight_was_updated() );
        //CPPUNIT_ASSERT_EQUAL(1.0, new_state.weight());
        //CPPUNIT_ASSERT(new_state.previous_state == NULL);
        
        //delete new_state.writable_model();
    }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestForestState );

