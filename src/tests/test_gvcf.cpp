#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/variantReader.hpp"
#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestGVCF : public CppUnit::TestCase {
  
    CPPUNIT_TEST_SUITE( TestGVCF );
    
    CPPUNIT_TEST( test_Constructors ); 
    CPPUNIT_TEST( test_empty_file ); 
    //CPPUNIT_TEST( test_new_block ); 
    //CPPUNIT_TEST( test_empty_block );
    //CPPUNIT_TEST( test_reset_data_to_first_entry );
    //CPPUNIT_TEST( test_read_NA12878 );    
    CPPUNIT_TEST_SUITE_END();

 private:
    VariantReader *gvcf_file;
    
 public:
    
    void test_empty_file(){
        gvcf_file = new VariantReader( "", GVCF );
        int testing_even_interval = 1000 ;
        CPPUNIT_ASSERT_NO_THROW( gvcf_file->even_interval_ = testing_even_interval );
        CPPUNIT_ASSERT_EQUAL( string(""), gvcf_file->file_name_ );
        CPPUNIT_ASSERT_EQUAL( size_t(-1), gvcf_file->file_length_ );  // As this is an empty file, it will reach the end of the file, and seekg will go to "infinity", which is size_t(-1)

        CPPUNIT_ASSERT_EQUAL( size_t(0), gvcf_file->nsam() );
        CPPUNIT_ASSERT_EQUAL( size_t(0), gvcf_file->nfield() );
        CPPUNIT_ASSERT_EQUAL( int(-1), gvcf_file->chrom() );
        CPPUNIT_ASSERT_EQUAL( EMPTY, gvcf_file->FileType );
        
        CPPUNIT_ASSERT_EQUAL( true, gvcf_file->eof_ ); // Empty file directly goes to the end of the file, as it may extend virtual data, it is not the end of the data
        CPPUNIT_ASSERT_EQUAL( false, gvcf_file->end_data() ); // Empty file directly goes to the end of the file, as it may extend virtual data, it is not the end of the data
        
        CPPUNIT_ASSERT_EQUAL( size_t(0), gvcf_file->current_line_index_ );
        
        CPPUNIT_ASSERT_EQUAL( size_t(0), gvcf_file->empty_file_line_counter_ );        
        CPPUNIT_ASSERT_EQUAL( int(0), gvcf_file->site() );
        
        // Start reading new "lines" in empty file
        for (size_t i = 0; i < 11; i++){
            CPPUNIT_ASSERT_NO_THROW( gvcf_file->read_new_line() );
            //CPPUNIT_ASSERT_EQUAL( size_t(0), gvcf_file->current_block_line_ ); // No need to initialize current_block_line_
            CPPUNIT_ASSERT_EQUAL( size_t(1+i), gvcf_file->empty_file_line_counter_ );
            CPPUNIT_ASSERT_EQUAL( int(i*testing_even_interval), gvcf_file->site() );
            }
            
        CPPUNIT_ASSERT_EQUAL( true, gvcf_file->end_data() );
        
        delete gvcf_file;
        }  
    
    void test_Constructors() {
        gvcf_file = new VariantReader("sim-1Samples2msdata1.gvcf", GVCF);
        CPPUNIT_ASSERT_EQUAL(string("sim-1Samples2msdata1.gvcf"), gvcf_file->file_name_);
        CPPUNIT_ASSERT_EQUAL(size_t(784647), gvcf_file->file_length_);  // wc test2sample.vcf, length of the file is 3227 characters

        CPPUNIT_ASSERT_EQUAL(size_t(1), gvcf_file->nsam()); // 1 diploid sample
        CPPUNIT_ASSERT_EQUAL(size_t(10), gvcf_file->nfield());
        CPPUNIT_ASSERT_EQUAL( GVCF, gvcf_file->FileType );
        
        CPPUNIT_ASSERT_EQUAL((int)-1, gvcf_file->site());
        CPPUNIT_ASSERT_EQUAL((int)-1, gvcf_file->chrom());
        
        
        CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
         *1	1	.	.	.	0	REFCALL;	END=381;	GT	./.
         */
        //CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(2), gvcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[0]);
        CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[1]);
        
        CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
         *1	382	rs0	A	T	67	PASS	NS=2;	GT	0|1
         */
        CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(int(382), gvcf_file->site());        
        CPPUNIT_ASSERT_EQUAL(size_t(2), gvcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[0]);
        CPPUNIT_ASSERT(gvcf_file->vec_of_sample_alt_bool[1]);
        
        delete gvcf_file;
        
        
        gvcf_file = new VariantReader("tests/test4sample.vcf", VCF);
        CPPUNIT_ASSERT_EQUAL(string("tests/test4sample.vcf"), gvcf_file->file_name_);
        CPPUNIT_ASSERT_EQUAL(size_t(6518), gvcf_file->file_length_);  // wc test2sample.vcf, length of the file is 6518 characters

        CPPUNIT_ASSERT_EQUAL(size_t(2), gvcf_file->nsam());
        CPPUNIT_ASSERT_EQUAL(size_t(11), gvcf_file->nfield());
        CPPUNIT_ASSERT_EQUAL( VCF, gvcf_file->FileType );
        
        CPPUNIT_ASSERT_EQUAL(int(-1), gvcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(-1), gvcf_file->chrom());
        cout<<"start here"<<endl;
cout<<gvcf_file->tmp_line<<endl;        
        CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2
         * 1	0	rs0	A	T	67	PASS	NS=2;	GT	0|0	1|1
         */
cout<<gvcf_file->tmp_line<<endl;
        CPPUNIT_ASSERT_EQUAL(int(0), gvcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(4), gvcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[0]);
        CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[1]);
        CPPUNIT_ASSERT(gvcf_file->vec_of_sample_alt_bool[2]);
        CPPUNIT_ASSERT(gvcf_file->vec_of_sample_alt_bool[3]);
        CPPUNIT_ASSERT_EQUAL(SEQ_INVARIANT, gvcf_file->prior_seq_state);
        CPPUNIT_ASSERT_EQUAL(SNP, gvcf_file->current_variant_state);
        
        CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2
         * 1       1       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|1
         */
        CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(4), gvcf_file->vec_of_sample_alt_bool.size());
        CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[0]);
        CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[1]);
        CPPUNIT_ASSERT(gvcf_file->vec_of_sample_alt_bool[2]);
        CPPUNIT_ASSERT(gvcf_file->vec_of_sample_alt_bool[3]); 
        
        delete gvcf_file;

        //gvcf_file = new VariantReader("tests/test6sample.vcf", GVCF);
        //CPPUNIT_ASSERT_EQUAL(string("tests/test6sample.vcf"), gvcf_file->file_name_);
        //CPPUNIT_ASSERT_EQUAL(size_t(13813), gvcf_file->file_length_);  // wc test2sample.vcf, length of the file is 13813 characters

        //CPPUNIT_ASSERT_EQUAL(size_t(3), gvcf_file->nsam());
        //CPPUNIT_ASSERT_EQUAL(size_t(12), gvcf_file->nfield());
        //CPPUNIT_ASSERT_EQUAL(SEQ_INVARIANT, gvcf_file->prior_seq_state);
        //CPPUNIT_ASSERT_EQUAL(SNP, gvcf_file->current_variant_state);
        
        //CPPUNIT_ASSERT_EQUAL(int(0), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(int(0), gvcf_file->chrom());
        
         ////cout << gvcf_file->buffer_lines[gvcf_file->current_block_line_]<<endl;   
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        ///*! read the following line
         //* #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2
         //* 1	0	rs0	A	T	67	PASS	NS=2;	GT	0|0	1|0	0|0
         //*/
        //CPPUNIT_ASSERT_EQUAL(int(0), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        //CPPUNIT_ASSERT_EQUAL(size_t(6), gvcf_file->vec_of_sample_alt_bool.size());
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[0]); 
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[1]);
        //CPPUNIT_ASSERT(gvcf_file->vec_of_sample_alt_bool[2]);
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[3]);   
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[4]); 
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[5]);      
        //CPPUNIT_ASSERT_EQUAL(SEQ_INVARIANT, gvcf_file->prior_seq_state);
        //CPPUNIT_ASSERT_EQUAL(SNP, gvcf_file->current_variant_state);

        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
                
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        ///*! read the following line
         //* #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2
         //* 1       4       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
         //*/
        //CPPUNIT_ASSERT_EQUAL(int(4), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        //CPPUNIT_ASSERT_EQUAL(size_t(6), gvcf_file->vec_of_sample_alt_bool.size());
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[0]);
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[1]);
        //CPPUNIT_ASSERT(gvcf_file->vec_of_sample_alt_bool[2]);
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[3]);
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[4]);
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[5]); 
        
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        ///*! read the following line
         //* #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2
         //*  1       7       rs0     A       T       67      PASS    NS=2;   GT      0|0     0|1     0|0
         //*/
        //CPPUNIT_ASSERT_EQUAL(int(7), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        //CPPUNIT_ASSERT_EQUAL(size_t(6), gvcf_file->vec_of_sample_alt_bool.size());
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[0]);
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[1]);
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[2]);
        //CPPUNIT_ASSERT(gvcf_file->vec_of_sample_alt_bool[3]); 
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[4]);
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[5]);
        
        //delete gvcf_file;

    }
    
    //void test_empty_block(){
        //int buff_len = 10;
        //gvcf_file = new VariantReader("tests/test6sample.vcf", GVCF, buff_len);
        //CPPUNIT_ASSERT_EQUAL(string("tests/test6sample.vcf"), gvcf_file->file_name_);
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line()); // vcf data line 1
        //CPPUNIT_ASSERT_EQUAL(size_t(buff_len), gvcf_file->buffer_lines.size());
        
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->empty_block()); 
        //CPPUNIT_ASSERT_EQUAL(size_t(0), gvcf_file->buffer_lines.size());
        
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_block()); 
        //CPPUNIT_ASSERT_EQUAL(size_t(buff_len), gvcf_file->buffer_lines.size());
        
        //delete gvcf_file;
    //}
    
    //void test_new_block(){
        //int buff_len = 3;
        //gvcf_file = new VariantReader("tests/test6sample.vcf", GVCF, buff_len);
        //CPPUNIT_ASSERT_EQUAL(string("tests/test6sample.vcf"), gvcf_file->file_name_);
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line()); // vcf data line 1
////1       0       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
////1       0       rs0     A       T       67      PASS    NS=2;   GT      0|1     0|0     1|0
////1       0       rs0     A       T       67      PASS    NS=2;   GT      0|1     0|0     1|0
////1       1       rs0     A       T       67      PASS    NS=2;   GT      0|0     0|0     1|0 *********
        //CPPUNIT_ASSERT_EQUAL(int(0), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(size_t(1), gvcf_file->current_block_line_); // line [0] is good, current_block_line is 0+1
        
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line()); 
////1       1       rs0     A       T       67      PASS    NS=2;   GT      0|0     0|0     1|0(repeat of the last line of the previous block)
////1       1       rs0     A       T       67      PASS    NS=2;   GT      0|1     0|0     1|0
////1       1       rs0     A       T       67      PASS    NS=2;   GT      0|1     0|0     1|0
////1       2       rs0     A       T       67      PASS    NS=2;   GT      1|0     0|1     0|1 *********
        //CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(size_t(1), gvcf_file->current_block_line_); // line [0] is good, current_block_line is 0+1

        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line()); 
////1       2       rs0     A       T       67      PASS    NS=2;   GT      1|0     0|1     0|1(repeat of the last line of the previous block)
////1       3       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
////1       3       rs0     A       T       67      PASS    NS=2;   GT      0|1     0|0     1|0
////1       3       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0 *********        
        //CPPUNIT_ASSERT_EQUAL(int(2), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(size_t(1), gvcf_file->current_block_line_); // line [0] is good, current_block_line is 0+1

        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());         
        //CPPUNIT_ASSERT_EQUAL(int(3), gvcf_file->site()); 
        //CPPUNIT_ASSERT_EQUAL(size_t(2), gvcf_file->current_block_line_); // line [1] is good, current_block_line is 1+1

        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line()); 
////1       3       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0(repeat of the last line of the previous block)
////1       4       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
////1       5       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
////1       5       rs0     A       T       67      PASS    NS=2;   GT      1|1     0|1     1|1 ********
        //CPPUNIT_ASSERT_EQUAL(int(4), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(size_t(2), gvcf_file->current_block_line_); // line [1] is good, current_block_line is 1+1

        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line()); 
        //CPPUNIT_ASSERT_EQUAL(int(5), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(size_t(0), gvcf_file->current_block_line_); // line [2] is good, current_block_line is 2+1, which is one less the size of the block, read a new block

////1       5       rs0     A       T       67      PASS    NS=2;   GT      1|1     0|1     1|1(repeat of the last line of the previous block)
////1       5       rs0     A       T       67      PASS    NS=2;   GT      1|1     0|1     1|1
////1       5       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
////1       5       rs0     A       T       67      PASS    NS=2;   GT      1|1     0|1     1|1 ********
////1       5       rs0     A       T       67      PASS    NS=2;   GT      1|1     0|1     1|1(repeat of the last line of the previous block)
////1       5       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
////1       5       rs0     A       T       67      PASS    NS=2;   GT      1|1     0|1     1|1
////1       5       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0 ********

        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line()); // vcf data line 7
////1       5       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0(repeat of the last line of the previous block)
////1       5       rs0     A       T       67      PASS    NS=2;   GT      0|0     0|1     0|0
////1       5       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
////1       5       rs0     A       T       67      PASS    NS=2;   GT      1|0     0|0     0|1 ********
////1       5       rs0     A       T       67      PASS    NS=2;   GT      1|0     0|0     0|1(repeat of the last line of the previous block)
////1       5       rs0     A       T       67      PASS    NS=2;   GT      1|0     0|0     0|1
////1       5       rs0     A       T       67      PASS    NS=2;   GT      0|0     0|1     0|0
////1       6       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0 ********

////1       6       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0(repeat of the last line of the previous block)
////1       6       rs0     A       T       67      PASS    NS=2;   GT      1|1     0|1     1|1
////1       6       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
////1       6       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0 ********

        //CPPUNIT_ASSERT_EQUAL(int(6), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(size_t(1), gvcf_file->current_block_line_); // line [0] is good, current_block_line is 1+1
////1       6       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0(repeat of the last line of the previous block)
////1       6       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
////1       6       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0
////1       6       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0 ********
////1       6       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0(repeat of the last line of the previous block)
////1       6       rs0     A       T       67      PASS    NS=2;   GT      1|1     0|1     1|1
////1       6       rs0     A       T       67      PASS    NS=2;   GT      1|1     0|1     1|1
////1       7       rs0     A       T       67      PASS    NS=2;   GT      0|0     0|1     0|0 ********

////1       7       rs0     A       T       67      PASS    NS=2;   GT      0|0     1|0     0|0(repeat of the last line of the previous block)
////1       7       rs0     A       T       67      PASS    NS=2;   GT      0|0     0|0     1|0
////1       8       rs0     A       T       67      PASS    NS=2;   GT      1|1     0|0     0|1
////1       8       rs0     A       T       67      PASS    NS=2;   GT      0|1     0|0     0|0 ********
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line()); // vcf data line 8
        //CPPUNIT_ASSERT_EQUAL(int(7), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(size_t(1), gvcf_file->current_block_line_); // line [0] is good, current_block_line is 1+1
        
        //delete gvcf_file;
    //}
    

    
    //void test_reset_data_to_first_entry(){
        //gvcf_file = new VariantReader("tests/test6sample.vcf", VCF);
        //CPPUNIT_ASSERT_EQUAL(string("tests/test6sample.vcf"), gvcf_file->file_name_);
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        ///*! read the following line
         //* #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2
         //* 1	0	rs0	A	T	67	PASS	NS=2;	GT	0|0	1|0	0|0
         //*/
        //CPPUNIT_ASSERT_EQUAL(int(0), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        //CPPUNIT_ASSERT_EQUAL(size_t(6), gvcf_file->vec_of_sample_alt_bool.size());
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[0]);
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[1]);
        //CPPUNIT_ASSERT(gvcf_file->vec_of_sample_alt_bool[2]);
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[3]);   
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[4]); 
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[5]);      
        
        //do{ 
            //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
            ////cout<<"Valid site "<<gvcf_file->site()<<endl;
        //}while(!gvcf_file->end_data());
        
        
        
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->reset_data_to_first_entry());
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        ///*! read the following line
         //* #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2
         //* 1	0	rs0	A	T	67	PASS	NS=2;	GT	0|0	1|0	0|0
         //*/
        //CPPUNIT_ASSERT_EQUAL(int(0), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        //CPPUNIT_ASSERT_EQUAL(size_t(6), gvcf_file->vec_of_sample_alt_bool.size());
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[0]);
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[1]);
        //CPPUNIT_ASSERT(gvcf_file->vec_of_sample_alt_bool[2]);
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[3]);   
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[4]); 
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[5]);      
                
        //delete gvcf_file;
        //}  
    
    //void test_read_NA18501(){
    //// CHECK if file is available 
    //ifstream my_file("tests/NA18501_CHROM1.vcf");
    //if ( !my_file.good() ){
        //cout << endl << endl << "tests/NA18501_CHROM1.vcf does not exist, unit test will fail" << endl << endl;
        //}
        //gvcf_file = new VariantReader("tests/NA18501_CHROM1.vcf", VCF);
        //gvcf_file->filter_window_ = 100;
        //CPPUNIT_ASSERT_EQUAL(string("tests/NA18501_CHROM1.vcf"), gvcf_file->file_name_);
        //// Just initialized, read buffer line: 
        ////1	52238	rs150021059	T	G	100	PASS	.	GT	1|0
        //// but not yet processed        
        //CPPUNIT_ASSERT_EQUAL(int(0), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(int(0), gvcf_file->chrom());
        //CPPUNIT_ASSERT_EQUAL(size_t(0), gvcf_file->vec_of_sample_alt_bool.size());
        
        ////cout << gvcf_file->buffer_lines[gvcf_file->current_block_line_]<<endl;        
        ////extract from line 
        ////1	52238	rs150021059	T	G	100	PASS	.	GT	1|0
        //CPPUNIT_ASSERT_NO_THROW( gvcf_file->read_new_line() );
        //CPPUNIT_ASSERT_EQUAL(int(52238), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        //CPPUNIT_ASSERT_EQUAL(size_t(2), gvcf_file->vec_of_sample_alt_bool.size());
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[1]);
        //CPPUNIT_ASSERT(gvcf_file->vec_of_sample_alt_bool[0]);
        ////cout << gvcf_file->buffer_lines[gvcf_file->current_block_line_]<<endl;
        ////extract from line 
        ////1	55164	rs3091274	C	A	100	PASS	.	GT	1|0
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        //CPPUNIT_ASSERT_EQUAL(int(55164), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        //CPPUNIT_ASSERT_EQUAL(size_t(2), gvcf_file->vec_of_sample_alt_bool.size());
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[1]);
        //CPPUNIT_ASSERT(gvcf_file->vec_of_sample_alt_bool[0]);

        ////cout << gvcf_file->buffer_lines[gvcf_file->current_block_line_]<<endl;
        ////extract from line 
        ////1	55299	rs10399749	C	T	100	PASS	.	GT	1|0
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        //CPPUNIT_ASSERT_EQUAL(int(55299), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        //CPPUNIT_ASSERT_EQUAL(size_t(2), gvcf_file->vec_of_sample_alt_bool.size());
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[1]);
        //CPPUNIT_ASSERT(gvcf_file->vec_of_sample_alt_bool[0]);


        ////cout << gvcf_file->buffer_lines[gvcf_file->current_block_line_]<<endl;
        ////extract from line 
        ////1	57952	rs189727433	A	C	100	PASS	.	GT	1|0
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        //CPPUNIT_ASSERT_EQUAL(int(57952), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        //CPPUNIT_ASSERT_EQUAL(size_t(2), gvcf_file->vec_of_sample_alt_bool.size());
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[1]);
        //CPPUNIT_ASSERT(gvcf_file->vec_of_sample_alt_bool[0]);

        ////cout << gvcf_file->buffer_lines[gvcf_file->current_block_line_]<<endl;
        ////extract from line 
        ////1	58814	rs114420996	G	A	100	PASS	.	GT	1|0
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        //CPPUNIT_ASSERT_EQUAL(int(58814), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        //CPPUNIT_ASSERT_EQUAL(size_t(2), gvcf_file->vec_of_sample_alt_bool.size());
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[1]);
        //CPPUNIT_ASSERT(gvcf_file->vec_of_sample_alt_bool[0]);

        ////cout << gvcf_file->buffer_lines[gvcf_file->current_block_line_]<<endl;
        ////extract from line 
        ////1	63671	rs116440577	G	A	100	PASS	.	GT	1|0
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        //CPPUNIT_ASSERT_EQUAL(int(63671), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        //CPPUNIT_ASSERT_EQUAL(size_t(2), gvcf_file->vec_of_sample_alt_bool.size());
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[1]);
        //CPPUNIT_ASSERT(gvcf_file->vec_of_sample_alt_bool[0]);


////cout << gvcf_file->buffer_lines[gvcf_file->current_block_line_]<<endl;
        ////Skip line
        ////1	63735	rs201888535	CCTA	C	100	PASS	.	GT	0|1
        ////extract from line 
        ////1	66176	rs28552463	T	A	100	PASS	.	GT	1|0
        //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        //CPPUNIT_ASSERT_EQUAL(int(66176), gvcf_file->site());
        //CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        //CPPUNIT_ASSERT_EQUAL(size_t(2), gvcf_file->vec_of_sample_alt_bool.size());
        //CPPUNIT_ASSERT(!gvcf_file->vec_of_sample_alt_bool[1]);
        //CPPUNIT_ASSERT(gvcf_file->vec_of_sample_alt_bool[0]);


        //gvcf_file->missing_data_threshold_ = 1000000;
        //do{             

            //CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
            
            //if ( gvcf_file->site() == 142565398 ){
                //CPPUNIT_ASSERT_EQUAL(MISSING, gvcf_file->prior_seq_state);
                //CPPUNIT_ASSERT_EQUAL(SNP, gvcf_file->current_variant_state);
                //}
            //else {
                //CPPUNIT_ASSERT_EQUAL(SEQ_INVARIANT, gvcf_file->prior_seq_state);
                ////CPPUNIT_ASSERT_EQUAL(SNP, gvcf_file->current_variant_state); // This not always true now! it can handle other polymophisms....
                //}
            
            //if (gvcf_file->site() > 142872424){
                //break;}
        //}while(!gvcf_file->end_data());
        
        //delete gvcf_file;
        //}
    
    
    
    };

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestGVCF );
