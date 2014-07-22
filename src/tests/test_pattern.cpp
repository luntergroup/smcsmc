#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "../src/pattern.hpp"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestPattern : public CppUnit::TestCase {
  
    CPPUNIT_TEST_SUITE( TestPattern );
    
//    CPPUNIT_TEST( test_extract_NumberOfSegment ); 
    
    CPPUNIT_TEST( test_for_faulty_pattern );

    CPPUNIT_TEST_SUITE_END();

    private:
        Pattern * pattern_tester;
        //Vcf *vcf_file;
        string pattern;
        double top_t ;
        //vector <size_t> seg_level1_vec_;
        //vector <size_t> seg_level2_vec_;
        vector <double> t_i;
        vector <double> new_ti;
        //size_t num_seg ;
        //const char * expr;
    
    public:
        void test_for_faulty_pattern () {
            top_t = 2.0; 
            
            pattern = "1*";
            
            pattern = "1adf*";
            
            pattern = "10-";
            pattern = "-";
            pattern = "+";
            pattern = "-5";
            pattern = "+4";
            
            }
    
    
    
        void test_extract_NumberOfSegment( ){
            top_t = 2.0; 
            
            pattern = "";
            
            pattern = "0";
            
            pattern = "1";
            
            pattern = "2";
            
            pattern = "3";
            
            pattern = "1+1";
            
            pattern = "3+2+2+2+2+2+3";
            //expr = pattern.c_str();
            pattern_tester = new Pattern (pattern, top_t);
            //num_seg = extract_NumberOfSegment( expr , seg_level1_vec_, seg_level2_vec_) ; 
            CPPUNIT_ASSERT_EQUAL( (size_t)16, pattern_tester->extract_NumberOfSegment() );
            t_i = extract_Segment( num_seg, this->top_t);
            CPPUNIT_ASSERT_EQUAL(top_t, t_i.back());
            CPPUNIT_ASSERT_EQUAL(t_i.size(), num_seg);            
            new_ti = regroup_Segment (t_i, seg_level1_vec_, seg_level2_vec_);
            CPPUNIT_ASSERT_EQUAL((size_t)7, new_ti.size());
            
            //seg_level1_vec_.clear();
            //seg_level2_vec_.clear();
            
            pattern = "3+2+2+2+2+2+2";
            expr = pattern.c_str();
            num_seg = extract_NumberOfSegment( expr , seg_level1_vec_, seg_level2_vec_) ; 
            CPPUNIT_ASSERT_EQUAL((size_t)15, num_seg);
            t_i = extract_Segment( num_seg, top_t);
            CPPUNIT_ASSERT_EQUAL(top_t, t_i.back());
            CPPUNIT_ASSERT_EQUAL(t_i.size(), num_seg);
            new_ti = regroup_Segment (t_i, seg_level1_vec_, seg_level2_vec_);
            CPPUNIT_ASSERT_EQUAL((size_t)7, new_ti.size());
            //seg_level1_vec_.clear();
            //seg_level2_vec_.clear();
            
            pattern = "1*4+25*2+1*4+1*6";
            expr = pattern.c_str();
            num_seg = extract_NumberOfSegment( expr , seg_level1_vec_, seg_level2_vec_) ; 
            CPPUNIT_ASSERT_EQUAL((size_t)64, num_seg);
            t_i = extract_Segment( num_seg, top_t);
            CPPUNIT_ASSERT_EQUAL(top_t, t_i.back());
            CPPUNIT_ASSERT_EQUAL(t_i.size(), num_seg);
            new_ti = regroup_Segment (t_i, seg_level1_vec_, seg_level2_vec_);
            CPPUNIT_ASSERT_EQUAL((size_t)28, new_ti.size());
            //seg_level1_vec_.clear();
            //seg_level2_vec_.clear();
            
            
            }  
    };
    
    
//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestPattern );
