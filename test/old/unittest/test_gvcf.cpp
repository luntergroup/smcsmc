#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "variantReader.hpp"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestGVCF : public CppUnit::TestCase {
  
    CPPUNIT_TEST_SUITE( TestGVCF );
    
    CPPUNIT_TEST( test_Constructors ); 
    CPPUNIT_TEST( test_empty_file ); 
    CPPUNIT_TEST( test_new_block ); 
    CPPUNIT_TEST( test_empty_block );
    CPPUNIT_TEST( test_reset_data_to_first_entry );
    CPPUNIT_TEST_SUITE_END();

 private:
    VariantReader *gvcf_file;
    
 public:
    
    void test_empty_file(){
        //VariantReader( "random", GVCF );
        CPPUNIT_ASSERT_THROW ( VariantReader( "random", GVCF ) , std::invalid_argument );
        //delete gvcf_file;

        gvcf_file = new VariantReader( "", GVCF );
        CPPUNIT_ASSERT_EQUAL( EMPTY, gvcf_file->FileType );
        int testing_even_interval = 1000 ;
        CPPUNIT_ASSERT_NO_THROW( gvcf_file->even_interval_ = testing_even_interval );
        CPPUNIT_ASSERT_EQUAL( string(""), gvcf_file->file_name_ );
        CPPUNIT_ASSERT_EQUAL( size_t(0), gvcf_file->file_length_ );  

        CPPUNIT_ASSERT_EQUAL( size_t(0), gvcf_file->nsam() );
        //CPPUNIT_ASSERT_EQUAL( size_t(0), gvcf_file->nfield() );
        CPPUNIT_ASSERT_EQUAL((int)0, gvcf_file->site());
        CPPUNIT_ASSERT_EQUAL((int)0, gvcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL((int)0, gvcf_file->previous_site_at_);
        CPPUNIT_ASSERT_EQUAL((int)0, gvcf_file->pervious_chrom_);

        
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
            //CPPUNIT_ASSERT_EQUAL( int(i*testing_even_interval), gvcf_file->site() );
            }
            
        CPPUNIT_ASSERT_EQUAL( true, gvcf_file->end_data() );
        
        delete gvcf_file;
        }  
    
    void test_Constructors() {
        gvcf_file = new VariantReader("sim-1Samples2msdata1.gvcf", GVCF);
        CPPUNIT_ASSERT_EQUAL(string("sim-1Samples2msdata1.gvcf"), gvcf_file->file_name_);
        CPPUNIT_ASSERT_EQUAL(size_t(784263), gvcf_file->file_length_);  // wc sim-1Samples2msdata1.vcf, length of the file is 784647 characters

        CPPUNIT_ASSERT_EQUAL(size_t(1), gvcf_file->nsam()); // 1 diploid sample
        //CPPUNIT_ASSERT_EQUAL(size_t(10), gvcf_file->nfield());
        CPPUNIT_ASSERT_EQUAL( GVCF, gvcf_file->FileType );
        
        CPPUNIT_ASSERT_EQUAL((int)0, gvcf_file->site());
        CPPUNIT_ASSERT_EQUAL((int)0, gvcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL((int)0, gvcf_file->previous_site_at_);
        CPPUNIT_ASSERT_EQUAL((int)0, gvcf_file->pervious_chrom_);
        
        
        CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
         *1	1	.	.	.	0	REFCALL;	END=381;	GT	./.
         */
        CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(381), gvcf_file->seg_end_site());
        CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(2), gvcf_file->int_vec_of_sample_alt.size());
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[0]);
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[1]);
        
        CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
         *1	382	rs0	A	T	67	PASS	NS=2;	GT	0|1
         */
        CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        //cout<< gvcf_file->tmp_line<<endl;
        CPPUNIT_ASSERT_EQUAL(int(382), gvcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(382), gvcf_file->seg_end_site());
        CPPUNIT_ASSERT_EQUAL(size_t(2), gvcf_file->int_vec_of_sample_alt.size());
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[0]);
        CPPUNIT_ASSERT(gvcf_file->int_vec_of_sample_alt[1]);


        CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
         *1       383     .       .       .       0       REFCALL;        END=1906;       GT      ./.
         */
        CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        //cout<< gvcf_file->tmp_line<<endl;
        CPPUNIT_ASSERT_EQUAL(int(383), gvcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1906), gvcf_file->seg_end_site());
        CPPUNIT_ASSERT_EQUAL(size_t(2), gvcf_file->int_vec_of_sample_alt.size());
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[0]);
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[1]);
                
        delete gvcf_file;
        
        
        gvcf_file = new VariantReader("sim-1Samples4msdata1_edited.gvcf", GVCF);
        CPPUNIT_ASSERT_EQUAL(string("sim-1Samples4msdata1_edited.gvcf"), gvcf_file->file_name_);

        CPPUNIT_ASSERT_EQUAL(size_t(1611032), gvcf_file->file_length_);  // wc sim-1Samples4msdata1_first40_lines.vcf, length of the file is 784647 characters

        CPPUNIT_ASSERT_EQUAL(size_t(2), gvcf_file->nsam()); // 2 diploid sample
        //CPPUNIT_ASSERT_EQUAL(size_t(10), gvcf_file->nfield());
        CPPUNIT_ASSERT_EQUAL( GVCF, gvcf_file->FileType );
        
        CPPUNIT_ASSERT_EQUAL((int)0, gvcf_file->site());
        CPPUNIT_ASSERT_EQUAL((int)0, gvcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL((int)0, gvcf_file->previous_site_at_);
        CPPUNIT_ASSERT_EQUAL((int)0, gvcf_file->pervious_chrom_);
        
        
        CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
         *1	1	.	.	.	0	REFCALL;	END=12;	GT	./.	./.
         */
        CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(12), gvcf_file->seg_end_site());
        CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(4), gvcf_file->int_vec_of_sample_alt.size());
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[0]);
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[1]);
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[2]);
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[3]);
        
        CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
         *1	13	rs0	A	T	67	PASS	NS=2;	GT	0|1	1|1
         */
        CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(int(13), gvcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(13), gvcf_file->seg_end_site());
        CPPUNIT_ASSERT_EQUAL(size_t(4), gvcf_file->int_vec_of_sample_alt.size());
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[0]);
        CPPUNIT_ASSERT(gvcf_file->int_vec_of_sample_alt[1]);
        CPPUNIT_ASSERT(gvcf_file->int_vec_of_sample_alt[2]);
        CPPUNIT_ASSERT(gvcf_file->int_vec_of_sample_alt[3]);
        
        delete gvcf_file;
    }
    
    void test_empty_block(){
        int buff_len = 10;
        gvcf_file = new VariantReader("sim-1Samples4msdata1_edited.gvcf", GVCF, buff_len);
        CPPUNIT_ASSERT_EQUAL(string("sim-1Samples4msdata1_edited.gvcf"), gvcf_file->file_name_);
        CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line()); // vcf data line 1
        CPPUNIT_ASSERT_EQUAL(size_t(buff_len), gvcf_file->buffer_lines.size());
        
        CPPUNIT_ASSERT_NO_THROW(gvcf_file->empty_block()); 
        CPPUNIT_ASSERT_EQUAL(size_t(0), gvcf_file->buffer_lines.size());
        
        CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_block()); 
        CPPUNIT_ASSERT_EQUAL(size_t(buff_len), gvcf_file->buffer_lines.size());
        
        delete gvcf_file;
    }
    
    void test_new_block(){
        int buff_len = 3;
        gvcf_file = new VariantReader("sim-1Samples4msdata1_edited.gvcf", GVCF, buff_len);
        CPPUNIT_ASSERT_EQUAL(string("sim-1Samples4msdata1_edited.gvcf"), gvcf_file->file_name_);
        CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line()); // vcf data line 1
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
        
        delete gvcf_file;
    }
    

    
    void test_reset_data_to_first_entry(){
        gvcf_file = new VariantReader("sim-1Samples4msdata1_edited.gvcf", VCF);
        CPPUNIT_ASSERT_EQUAL(string("sim-1Samples4msdata1_edited.gvcf"), gvcf_file->file_name_);
        
        CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2
         * 1       1       .       .       .       0       REFCALL;        END=12; GT      ./.     ./.
         */
        CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(4), gvcf_file->int_vec_of_sample_alt.size());
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[0]);
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[1]);
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[2]);
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[3]);   
        
        
        do{ 
            CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
            //cout<<"Valid site "<<gvcf_file->site()<<endl;
        }while(!gvcf_file->end_data());
        
        CPPUNIT_ASSERT_NO_THROW(gvcf_file->reset_data_to_first_entry());
        CPPUNIT_ASSERT_NO_THROW(gvcf_file->read_new_line());
        /*! read the following line
         * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2
         * 1	0	rs0	A	T	67	PASS	NS=2;	GT	0|0	1|0	0|0
         */
        CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->site());
        CPPUNIT_ASSERT_EQUAL(int(1), gvcf_file->chrom());
        CPPUNIT_ASSERT_EQUAL(size_t(4), gvcf_file->int_vec_of_sample_alt.size());
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[0]);
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[1]);
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[2]);
        CPPUNIT_ASSERT(!gvcf_file->int_vec_of_sample_alt[3]);   
                
        delete gvcf_file;
        }  
        
    };

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestGVCF );
