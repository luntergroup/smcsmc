#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <boost/lexical_cast.hpp> 

#include "../src/vcf.hpp"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestVcf : public CppUnit::TestCase {
  
    CPPUNIT_TEST_SUITE( TestVcf );
    
    CPPUNIT_TEST( test_Constructors ); 
    CPPUNIT_TEST( test_empty_file ); 
    CPPUNIT_TEST( test_new_block ); 
    CPPUNIT_TEST( test_empty_block );
    CPPUNIT_TEST( test_reset_VCF_to_data );
    CPPUNIT_TEST( test_read_NA18501 );    
    CPPUNIT_TEST_SUITE_END();

 private:
    Vcf *vcf_file;
    
 public:
    
    void test_empty_file(){
        vcf_file = new Vcf( "" );
        double testing_even_interval = 1000 ;
        CPPUNIT_ASSERT_NO_THROW( vcf_file->even_interval_ = testing_even_interval );
        CPPUNIT_ASSERT_EQUAL( string(""), vcf_file->file_name_ );
        CPPUNIT_ASSERT_EQUAL( size_t(-1), vcf_file->vcf_length_ );  // As this is an empty file, it will reach the end of the file, and seekg will go to "infinity", which is size_t(-1)

        CPPUNIT_ASSERT_EQUAL( size_t(0), vcf_file->nsam() );
        CPPUNIT_ASSERT_EQUAL( size_t(0), vcf_file->nfield() );
        CPPUNIT_ASSERT_EQUAL( int(0), vcf_file->chrom() );
        CPPUNIT_ASSERT_EQUAL( false, vcf_file->withdata() );
        
        CPPUNIT_ASSERT_EQUAL( true, vcf_file->eof() ); // Empty file directly goes to the end of the file, as it may extend virtual data, it is not the end of the data
        CPPUNIT_ASSERT_EQUAL( false, vcf_file->end_data() ); // Empty file directly goes to the end of the file, as it may extend virtual data, it is not the end of the data
        
        CPPUNIT_ASSERT_EQUAL( size_t(0), vcf_file->current_line_index_ );
        
        CPPUNIT_ASSERT_EQUAL( size_t(0), vcf_file->empty_file_line_counter_ );        
        CPPUNIT_ASSERT_EQUAL( double(0), vcf_file->site() );
        
        // Start reading new "lines" in empty file
        for (size_t i = 0; i < 11; i++){
            CPPUNIT_ASSERT_NO_THROW( vcf_file->read_new_line() );
            CPPUNIT_ASSERT_EQUAL( size_t(0), vcf_file->current_block_line_ );
            CPPUNIT_ASSERT_EQUAL( size_t(1+i), vcf_file->empty_file_line_counter_ );
            CPPUNIT_ASSERT_EQUAL( double(i*testing_even_interval), vcf_file->site() );
            }
            
        CPPUNIT_ASSERT_EQUAL( true, vcf_file->end_data() );
        
        delete vcf_file;
        }  
    
    void test_Constructors() {
        vcf_file = new Vcf("tests/test2sample.vcf");
        CPPUNIT_ASSERT_EQUAL(string("tests/test2sample.vcf"), vcf_file->file_name_);
        CPPUNIT_ASSERT_EQUAL(size_t(3227), vcf_file->vcf_length_);  // wc test2sample.vcf, length of the file is 3227 characters

        CPPUNIT_ASSERT_EQUAL(size_t(1), vcf_file->nsam());
        CPPUNIT_ASSERT_EQUAL(size_t(10), vcf_file->nfield());
        CPPUNIT_ASSERT_EQUAL(true, vcf_file->withdata());
        
        CPPUNIT_ASSERT_EQUAL(double(0), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(0), vcf_file->chrom());
        
        
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
         * 1	1	rs0	A	T	67	PASS	NS=2;	GT	1|0
         */
        CPPUNIT_ASSERT_EQUAL(double(1), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), vcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(2), vcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[0]);
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[1]);
        
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
         * 1	2	rs0	A	T	67	PASS	NS=2;	GT	0|1
         */
        CPPUNIT_ASSERT_EQUAL(double(2), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), vcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(2), vcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[0]);
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[1]);
        
        delete vcf_file;
        
        
        vcf_file = new Vcf("tests/test4sample.vcf");
        CPPUNIT_ASSERT_EQUAL(string("tests/test4sample.vcf"), vcf_file->file_name_);
        CPPUNIT_ASSERT_EQUAL(size_t(6518), vcf_file->vcf_length_);  // wc test2sample.vcf, length of the file is 6518 characters

        CPPUNIT_ASSERT_EQUAL(size_t(2), vcf_file->nsam());
        CPPUNIT_ASSERT_EQUAL(size_t(11), vcf_file->nfield());
        CPPUNIT_ASSERT_EQUAL(true, vcf_file->withdata());
        
        CPPUNIT_ASSERT_EQUAL(double(0), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(0), vcf_file->chrom());
        
        
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2
         * 1	0	rs0	A	T	67	PASS	NS=2;	GT	0|0	1|1
         */
        CPPUNIT_ASSERT_EQUAL(double(0), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), vcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(4), vcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[0]);
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[1]);
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[2]);
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[3]);        
        
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2
         * 1       1       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|1
         */
        CPPUNIT_ASSERT_EQUAL(double(1), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), vcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(4), vcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[0]);
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[1]);
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[2]);
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[3]); 
        
        delete vcf_file;

        vcf_file = new Vcf("tests/test6sample.vcf");
        CPPUNIT_ASSERT_EQUAL(string("tests/test6sample.vcf"), vcf_file->file_name_);
        CPPUNIT_ASSERT_EQUAL(size_t(13813), vcf_file->vcf_length_);  // wc test2sample.vcf, length of the file is 13813 characters

        CPPUNIT_ASSERT_EQUAL(size_t(3), vcf_file->nsam());
        CPPUNIT_ASSERT_EQUAL(size_t(12), vcf_file->nfield());
        CPPUNIT_ASSERT_EQUAL(true, vcf_file->withdata());
        
        CPPUNIT_ASSERT_EQUAL(double(0), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(0), vcf_file->chrom());
        
         //cout << vcf_file->buffer_lines[vcf_file->current_block_line_]<<endl;   
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2
         * 1	0	rs0	A	T	67	PASS	NS=2;	GT	0|0	1|0	0|0
         */
        CPPUNIT_ASSERT_EQUAL(double(0), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), vcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(6), vcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[0]); 
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[1]);
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[2]);
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[3]);   
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[4]); 
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[5]);      
        
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
                
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2
         * 1       4       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
         */
        CPPUNIT_ASSERT_EQUAL(double(4), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), vcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(6), vcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[0]);
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[1]);
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[2]);
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[3]);
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[4]);
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[5]); 
        
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2
         *  1       7       rs0     A       T       67      PASS    NS=2;   GT      0|0     0|1     0|0
         */
        CPPUNIT_ASSERT_EQUAL(double(7), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), vcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(6), vcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[0]);
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[1]);
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[2]);
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[3]); 
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[4]);
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[5]);
        
        delete vcf_file;

    }
    
    void test_empty_block(){
        int buff_len = 10;
        vcf_file = new Vcf("tests/test6sample.vcf", buff_len);
        CPPUNIT_ASSERT_EQUAL(string("tests/test6sample.vcf"), vcf_file->file_name_);
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line()); // vcf data line 1
        CPPUNIT_ASSERT_EQUAL(size_t(buff_len), vcf_file->buffer_lines.size());
        
        CPPUNIT_ASSERT_NO_THROW(vcf_file->empty_block()); 
        CPPUNIT_ASSERT_EQUAL(size_t(0), vcf_file->buffer_lines.size());
        
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_block()); 
        CPPUNIT_ASSERT_EQUAL(size_t(buff_len), vcf_file->buffer_lines.size());
        
        delete vcf_file;
    }
    
    void test_new_block(){
        int buff_len = 3;
        vcf_file = new Vcf("tests/test6sample.vcf", buff_len);
        CPPUNIT_ASSERT_EQUAL(string("tests/test6sample.vcf"), vcf_file->file_name_);
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line()); // vcf data line 1
        //1       0       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
        //1       0       rs0     A       T       67      PASS    NS=2;   GT      0|1     0|0     1|0
        //1       0       rs0     A       T       67      PASS    NS=2;   GT      0|1     0|0     1|0
        CPPUNIT_ASSERT_EQUAL(double(0), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(size_t(1), vcf_file->current_block_line_);
        
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line()); 
//1       1       rs0     A       T       67      PASS    NS=2;   GT      0|0     0|0     1|0
//1       1       rs0     A       T       67      PASS    NS=2;   GT      0|1     0|0     1|0
//1       1       rs0     A       T       67      PASS    NS=2;   GT      0|1     0|0     1|0
        CPPUNIT_ASSERT_EQUAL(double(1), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(size_t(1), vcf_file->current_block_line_);

        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line()); 
//1       2       rs0     A       T       67      PASS    NS=2;   GT      1|0     0|1     0|1
//1       3       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
//1       3       rs0     A       T       67      PASS    NS=2;   GT      0|1     0|0     1|0
        CPPUNIT_ASSERT_EQUAL(double(2), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(size_t(1), vcf_file->current_block_line_); // new block has been read
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());         
        CPPUNIT_ASSERT_EQUAL(double(3), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(size_t(2), vcf_file->current_block_line_);

        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line()); 
//1       3       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
//1       4       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
//1       5       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
        CPPUNIT_ASSERT_EQUAL(double(4), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(size_t(2), vcf_file->current_block_line_);

        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line()); // vcf data line 6
//1       5       rs0     A       T       67      PASS    NS=2;   GT      1|1     0|1     1|1
//1       5       rs0     A       T       67      PASS    NS=2;   GT      1|1     0|1     1|1
//1       5       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
        CPPUNIT_ASSERT_EQUAL(double(5), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(size_t(0), vcf_file->current_block_line_); // new block has been read

        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line()); // vcf data line 7
//1       5       rs0     A       T       67      PASS    NS=2;   GT      1|1     0|1     1|1
//1       5       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
//1       5       rs0     A       T       67      PASS    NS=2;   GT      1|1     0|1     1|1
//1       5       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
//1       5       rs0     A       T       67      PASS    NS=2;   GT      0|0     0|1     0|0
//1       5       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
//1       5       rs0     A       T       67      PASS    NS=2;   GT      1|0     0|0     0|1
//1       5       rs0     A       T       67      PASS    NS=2;   GT      1|0     0|0     0|1
//1       5       rs0     A       T       67      PASS    NS=2;   GT      0|0     0|1     0|0

//1       6       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
//1       6       rs0     A       T       67      PASS    NS=2;   GT      1|1     0|1     1|1
//1       6       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
        CPPUNIT_ASSERT_EQUAL(double(6), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(size_t(1), vcf_file->current_block_line_);

//1       6       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
//1       6       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
//1       6       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
//1       6       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
//1       6       rs0     A       T       67      PASS    NS=2;   GT      1|1     0|1     1|1
//1       6       rs0     A       T       67      PASS    NS=2;   GT      1|1     0|1     1|1
// when reading in new lines, new blocks have been formed
//1       7       rs0     A       T       67      PASS    NS=2;   GT      0|0     0|1     0|0
//1       7       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
//1       7       rs0     A       T       67      PASS    NS=2;   GT      0|0     0|0     1|0
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line()); // vcf data line 8
        CPPUNIT_ASSERT_EQUAL(double(7), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(size_t(1), vcf_file->current_block_line_);
        
        delete vcf_file;
    }
    

    
    void test_reset_VCF_to_data(){
        vcf_file = new Vcf("tests/test6sample.vcf");
        CPPUNIT_ASSERT_EQUAL(string("tests/test6sample.vcf"), vcf_file->file_name_);
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2
         * 1	0	rs0	A	T	67	PASS	NS=2;	GT	0|0	1|0	0|0
         */
        CPPUNIT_ASSERT_EQUAL(double(0), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), vcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(6), vcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[0]);
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[1]);
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[2]);
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[3]);   
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[4]); 
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[5]);      
        
        do{ 
            CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
            //cout<<"Valid site "<<vcf_file->site()<<endl;
        }while(!vcf_file->end_data());
        
        
        
        CPPUNIT_ASSERT_NO_THROW(vcf_file->reset_VCF_to_data());
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2
         * 1	0	rs0	A	T	67	PASS	NS=2;	GT	0|0	1|0	0|0
         */
        CPPUNIT_ASSERT_EQUAL(double(0), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), vcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(6), vcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[0]);
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[1]);
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[2]);
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[3]);   
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[4]); 
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[5]);      
                
        delete vcf_file;
        }  
    
    void test_read_NA18501(){
        vcf_file = new Vcf("tests/NA18501_CHROM1.vcf");
        vcf_file->filter_window_ = 100;
        CPPUNIT_ASSERT_EQUAL(string("tests/NA18501_CHROM1.vcf"), vcf_file->file_name_);
        // Just initialized, read buffer line: 
        //1	52238	rs150021059	T	G	100	PASS	.	GT	1|0
        // but not yet processed        
        CPPUNIT_ASSERT_EQUAL(double(0), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(0), vcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(0), vcf_file->vec_of_sample_alt_bool.size());

        //cout << vcf_file->buffer_lines[vcf_file->current_block_line_]<<endl;        
        //extract from line 
        //1	52238	rs150021059	T	G	100	PASS	.	GT	1|0
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        CPPUNIT_ASSERT_EQUAL(double(52238), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), vcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(2), vcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[1]);
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[0]);
        
        //cout << vcf_file->buffer_lines[vcf_file->current_block_line_]<<endl;
        //extract from line 
        //1	55164	rs3091274	C	A	100	PASS	.	GT	1|0
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        CPPUNIT_ASSERT_EQUAL(double(55164), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), vcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(2), vcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[1]);
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[0]);

        //cout << vcf_file->buffer_lines[vcf_file->current_block_line_]<<endl;
        //extract from line 
        //1	55299	rs10399749	C	T	100	PASS	.	GT	1|0
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        CPPUNIT_ASSERT_EQUAL(double(55299), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), vcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(2), vcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[1]);
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[0]);


        //cout << vcf_file->buffer_lines[vcf_file->current_block_line_]<<endl;
        //extract from line 
        //1	57952	rs189727433	A	C	100	PASS	.	GT	1|0
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        CPPUNIT_ASSERT_EQUAL(double(57952), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), vcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(2), vcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[1]);
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[0]);

        //cout << vcf_file->buffer_lines[vcf_file->current_block_line_]<<endl;
        //extract from line 
        //1	58814	rs114420996	G	A	100	PASS	.	GT	1|0
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        CPPUNIT_ASSERT_EQUAL(double(58814), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), vcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(2), vcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[1]);
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[0]);

        //cout << vcf_file->buffer_lines[vcf_file->current_block_line_]<<endl;
        //extract from line 
        //1	63671	rs116440577	G	A	100	PASS	.	GT	1|0
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        CPPUNIT_ASSERT_EQUAL(double(63671), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), vcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(2), vcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[1]);
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[0]);


//cout << vcf_file->buffer_lines[vcf_file->current_block_line_]<<endl;
        //Skip line
        //1	63735	rs201888535	CCTA	C	100	PASS	.	GT	0|1
        //extract from line 
        //1	66176	rs28552463	T	A	100	PASS	.	GT	1|0
        CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
        CPPUNIT_ASSERT_EQUAL(double(66176), vcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), vcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(2), vcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!vcf_file->vec_of_sample_alt_bool[1]);
        CPPUNIT_ASSERT(vcf_file->vec_of_sample_alt_bool[0]);

        do{             

            CPPUNIT_ASSERT_NO_THROW(vcf_file->read_new_line());
            //cout<<"Valid site "<<vcf_file->site()<<endl;
            if (vcf_file->site() > 121426494){
                break;}
        }while(!vcf_file->end_data());
        
        delete vcf_file;
        }
    
    
    
    };

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestVcf );
