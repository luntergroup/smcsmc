#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "../src/pattern.hpp"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestPattern : public CppUnit::TestCase {
  
    CPPUNIT_TEST_SUITE( TestPattern );
    
    CPPUNIT_TEST( test_extract_NumberOfSegment ); 

    CPPUNIT_TEST_SUITE_END();

    private:
        //Vcf *vcf_file;
        string pattern;
        double top_t = 2.0;    
        vector <size_t> seg_level1_vec;
        vector <size_t> seg_level2_vec;
        vector <double> t_i;
        vector <double> new_ti;
        size_t num_seg ;
        const char * expr;
    
    public:
        void test_extract_NumberOfSegment( ){
            
            pattern = "3+2+2+2+2+2+3";
            expr = pattern.c_str();
            num_seg = extract_NumberOfSegment( expr , seg_level1_vec, seg_level2_vec) ; 
            CPPUNIT_ASSERT_EQUAL((size_t)16, num_seg);
            t_i = extract_Segment( num_seg, top_t);
            CPPUNIT_ASSERT_EQUAL(top_t, t_i.back());
            CPPUNIT_ASSERT_EQUAL(t_i.size(), num_seg);            
            new_ti = regroup_Segment (t_i, seg_level1_vec, seg_level2_vec);
            CPPUNIT_ASSERT_EQUAL((size_t)7, new_ti.size());
            
            seg_level1_vec.clear();
            seg_level2_vec.clear();
            
            pattern = "3+2+2+2+2+2+2";
            expr = pattern.c_str();
            num_seg = extract_NumberOfSegment( expr , seg_level1_vec, seg_level2_vec) ; 
            CPPUNIT_ASSERT_EQUAL((size_t)15, num_seg);
            t_i = extract_Segment( num_seg, top_t);
            CPPUNIT_ASSERT_EQUAL(top_t, t_i.back());
            CPPUNIT_ASSERT_EQUAL(t_i.size(), num_seg);
            new_ti = regroup_Segment (t_i, seg_level1_vec, seg_level2_vec);
            CPPUNIT_ASSERT_EQUAL((size_t)7, new_ti.size());
            seg_level1_vec.clear();
            seg_level2_vec.clear();
            
            pattern = "1*4+25*2+1*4+1*6";
            expr = pattern.c_str();
            num_seg = extract_NumberOfSegment( expr , seg_level1_vec, seg_level2_vec) ; 
            CPPUNIT_ASSERT_EQUAL((size_t)64, num_seg);
            t_i = extract_Segment( num_seg, top_t);
            CPPUNIT_ASSERT_EQUAL(top_t, t_i.back());
            CPPUNIT_ASSERT_EQUAL(t_i.size(), num_seg);
            new_ti = regroup_Segment (t_i, seg_level1_vec, seg_level2_vec);
            CPPUNIT_ASSERT_EQUAL((size_t)28, new_ti.size());
            seg_level1_vec.clear();
            seg_level2_vec.clear();
            
            
            }  
    };
    
    
//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestPattern );
