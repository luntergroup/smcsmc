#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <boost/lexical_cast.hpp> 

#include "count.hpp"
#include <algorithm>    // std::equal

#include "model.h"
//#include "param.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestCount : public CppUnit::TestCase {
  
    CPPUNIT_TEST_SUITE( TestCount );
    
    CPPUNIT_TEST( test_extract_change_time ); 
    CPPUNIT_TEST( test_setters_and_getters );
    CPPUNIT_TEST( test_Constructors ); 
    CPPUNIT_TEST_SUITE_END();

 private:
    Model * model_null ;
    CountModel * countNe_null;
    Model * model ;
    CountModel * countNe;

    double default_pop_size;
    void init(){
        model_null = new Model();    
        model = new Model();    
        char *argv[] = { "scrm", "2", "1", "-eN", "1.1", "2", "-eN", "3", "1", "-eN", "4.2", "2.3"};
        CPPUNIT_ASSERT_NO_THROW( Param(12, argv).parse(*model) );  
        CPPUNIT_ASSERT_NO_THROW(default_pop_size = model->default_pop_size);           
        
        //countNe_null = NULL;
        //countNe = NULL;
        
    }
    
    void finish(){
        if (model) delete model;
        if (model_null) delete model_null;    
        if (countNe_null) delete countNe_null;
        if (countNe) delete countNe;
    }

 public:

    void test_extract_change_time(){
        // 
        this->init();
        //model_null = new Model();    
        //model = new Model();    
        //char *argv[] = { "scrm", "2", "1", "-eN", "1.1", "2", "-eN", "3", "1", "-eN", "4.2", "2.3"};
        //CPPUNIT_ASSERT_NO_THROW( Param(12, argv).parse(*model) );  
        //CPPUNIT_ASSERT_NO_THROW(default_pop_size = model->default_pop_size);           
        
        vector<double> tested_time1 = extrcat_change_times_frm_model(model_null);
        double expected_time1[] = {0}; 
        CPPUNIT_ASSERT_EQUAL((size_t)1,tested_time1.size());
        CPPUNIT_ASSERT(std::equal(tested_time1.begin(), tested_time1.end(), expected_time1));


        vector<double> tested_time2 = extrcat_change_times_frm_model(model);
        
        CPPUNIT_ASSERT_EQUAL((size_t)4,tested_time2.size());
        double expected_time2[] = {0, 1.1 * 4 * default_pop_size, 3.0 * 4 * default_pop_size, 4.2 * 4 * default_pop_size}; 
        CPPUNIT_ASSERT(std::equal(tested_time2.begin(), tested_time2.end(), expected_time2));
        
        this->finish();
    }

    void test_Constructors() {
        this->init();
        countNe_null = new Count(model_null);
        countNe = new Count(model);        
        this->finish();
    }
    
    void test_setters_and_getters(){
        this->init();
        countNe_null = new Count(model_null);
        CPPUNIT_ASSERT_EQUAL((size_t)1, countNe_null->total_coal_count.size());
        CPPUNIT_ASSERT_EQUAL((size_t)1, countNe_null->total_weighted_coal_opportunity.size());
        CPPUNIT_ASSERT_EQUAL((size_t)1, countNe_null->Nehat.size());
        CPPUNIT_ASSERT_EQUAL((double)1/(2*default_pop_size), countNe_null->total_coal_count[0]);
        CPPUNIT_ASSERT_EQUAL((double)1, countNe_null->total_weighted_coal_opportunity[0]);
        CPPUNIT_ASSERT_EQUAL((double)default_pop_size, countNe_null->Nehat[0]);
        
        countNe = new Count(model); 
        CPPUNIT_ASSERT_EQUAL((size_t)4, countNe->total_coal_count.size());
        CPPUNIT_ASSERT_EQUAL((size_t)4, countNe->total_weighted_coal_opportunity.size());
        CPPUNIT_ASSERT_EQUAL((size_t)4, countNe->Nehat.size());
        
        CPPUNIT_ASSERT_EQUAL((double)1/(2*default_pop_size), countNe->total_coal_count[0]);
        CPPUNIT_ASSERT_EQUAL(1/(2*default_pop_size * 2), countNe->total_coal_count[1]);
        CPPUNIT_ASSERT_EQUAL(1/(2*default_pop_size * 1), countNe->total_coal_count[2]);
        CPPUNIT_ASSERT_EQUAL(1/(2*default_pop_size * 2.3), countNe->total_coal_count[3]);
        
        CPPUNIT_ASSERT_EQUAL((double)1, countNe->total_weighted_coal_opportunity[0]);
        CPPUNIT_ASSERT_EQUAL((double)1, countNe->total_weighted_coal_opportunity[1]);
        CPPUNIT_ASSERT_EQUAL((double)1, countNe->total_weighted_coal_opportunity[2]);
        CPPUNIT_ASSERT_EQUAL((double)1, countNe->total_weighted_coal_opportunity[3]);
        
        CPPUNIT_ASSERT_EQUAL(default_pop_size, countNe->Nehat[0]);
        CPPUNIT_ASSERT_EQUAL(default_pop_size*double(2), countNe->Nehat[1]);
        CPPUNIT_ASSERT_EQUAL(default_pop_size*double(1), countNe->Nehat[2]);
        CPPUNIT_ASSERT_EQUAL(default_pop_size*2.3,  countNe->Nehat[3]);
        
        this->finish();
    }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestCount );
