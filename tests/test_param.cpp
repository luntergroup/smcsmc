#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
//#include <boost/lexical_cast.hpp> 

#include "pfparam.hpp"
#include "help.hpp"
#include "variantReader.hpp"
#include "param.h"
#include "model.h"
#include "forest.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestParam : public CppUnit::TestCase {
  
  CPPUNIT_TEST_SUITE( TestParam );

  CPPUNIT_TEST( testParse );
  //CPPUNIT_TEST( testParseMigrationOptions );
  //CPPUNIT_TEST( testParseGrowthOptions );
  CPPUNIT_TEST( testReadInput );
  //CPPUNIT_TEST(test_initialize_model);
  
  
  CPPUNIT_TEST_SUITE_END();

 private:
    Model model;
    VariantReader* vcf_file;
    
 public:
  void testParse() {
    CPPUNIT_ASSERT_NO_THROW( Param().parse(model) );

    char *argv1[] = { "smcsmc"};
    PfParam pfARG_para1(1, argv1);
    CPPUNIT_ASSERT_EQUAL( (size_t)100, pfARG_para1.N );
    CPPUNIT_ASSERT_EQUAL( true, pfARG_para1.log_bool );  
    //CPPUNIT_ASSERT_EQUAL( 200, pfARG_para1.buff_length );  
    CPPUNIT_ASSERT_EQUAL( (double)0, pfARG_para1.lag );  
    CPPUNIT_ASSERT_EQUAL( string("pfARG"), pfARG_para1.out_NAME_prefix );  
    CPPUNIT_ASSERT_EQUAL( 0.5, pfARG_para1.ESS() );  
    CPPUNIT_ASSERT_EQUAL( true, pfARG_para1.ESS_default_bool );  
    CPPUNIT_ASSERT_EQUAL( false, pfARG_para1.online_bool );  
    CPPUNIT_ASSERT_EQUAL( (int)400, pfARG_para1.heat_seq_window );  
    CPPUNIT_ASSERT_EQUAL( 0, pfARG_para1.EM_steps );  
    CPPUNIT_ASSERT_EQUAL( false, pfARG_para1.EM_bool );  


    //Param pars = Param(pfARG_para1.scrm_argc_, pfARG_para1.scrm_argv_, false);    
    //CPPUNIT_ASSERT_EQUAL( false, pars.directly_called_);
    //CPPUNIT_ASSERT_NO_THROW( pars.parse(model) );
    ////CPPUNIT_ASSERT_EQUAL( false, pars.tree_bool );
    //CPPUNIT_ASSERT_EQUAL( false, pars.seg_bool() );  
    ////CPPUNIT_ASSERT_EQUAL( false, pars.tmrca_bool );


    char *argv2[] = { "smcsmc", "-Np", "10", "-log", "-lag", "10", "-online", "-EM", "5", "-ESS", "0.75"};
    PfParam pfARG_para2(11, argv2);
    CPPUNIT_ASSERT_EQUAL( (size_t)10, pfARG_para2.N );
    CPPUNIT_ASSERT_EQUAL( true, pfARG_para2.log_bool );  
    //CPPUNIT_ASSERT_EQUAL( (int)200, pfARG_para2.buff_length );  
    CPPUNIT_ASSERT_EQUAL( (double)10, pfARG_para2.lag );  
    CPPUNIT_ASSERT_EQUAL( (string)"pfARG", pfARG_para2.out_NAME_prefix );  
    
    CPPUNIT_ASSERT( pfARG_para2.ESS() == 0.75);  // DEBUG this line does not under valgrind ???
    //cout<<" !!!!!!!! "<<pfARG_para2.ESS<<endl;
    CPPUNIT_ASSERT_EQUAL( false, pfARG_para2.ESS_default_bool );  
    CPPUNIT_ASSERT_EQUAL( true, pfARG_para2.online_bool );  
    CPPUNIT_ASSERT_EQUAL( int(400), pfARG_para2.heat_seq_window );  
    CPPUNIT_ASSERT_EQUAL( 5, pfARG_para2.EM_steps );  
    CPPUNIT_ASSERT_EQUAL( true, pfARG_para2.EM_bool );  

    //Param pars2 = Param(7, pfARG_para2.scrm_argv_, false);    
    //CPPUNIT_ASSERT_EQUAL( false, pars2.directly_called_);
    //CPPUNIT_ASSERT_NO_THROW( pars2.parse(model) );
    ////CPPUNIT_ASSERT_EQUAL( false, pars2.tree_bool );
    //CPPUNIT_ASSERT_EQUAL( true, pars2.seg_bool() );  
    ////CPPUNIT_ASSERT_EQUAL( false, pars2.tmrca_bool );

    
    char *argv3[] = { "smcsmc", "-log", "-lag", "10", "-online" , "-Na"};
    //CPPUNIT_ASSERT_THROW( PfParam(6, argv3), std::invalid_argument); // "-Na" is invalid
    argv3[5] = "-ESS";
    CPPUNIT_ASSERT_THROW( PfParam(6, argv3), std::invalid_argument ); // "-ESS" expect more input
    argv3[5] = "-EM";
    CPPUNIT_ASSERT_THROW( PfParam(6, argv3), std::invalid_argument ); // "-EM" expect more input
    argv3[2] = "-EM";
    argv3[5] = "-lag";
    CPPUNIT_ASSERT_THROW( PfParam(6, argv3), std::invalid_argument ); // "-lag" expect more input
    argv3[5] = "-Np";
    CPPUNIT_ASSERT_THROW( PfParam(6, argv3), std::invalid_argument ); // "-Np" expect more input
    argv3[6] = "10.4";                                                      // "-Np" expect an integer
    //CPPUNIT_ASSERT_THROW( PfParam(7, argv3), boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::bad_lexical_cast> > ); 
    CPPUNIT_ASSERT_THROW( PfParam(7, argv3), std::invalid_argument ); // "-Np" expect more input
    argv3[6] = "10";
    CPPUNIT_ASSERT_NO_THROW( PfParam(7, argv3)); // "Valid input! "
    
    //char *argv4[] = { "smcsmc", "-logo"}; // "-logo" is unknown option
    //CPPUNIT_ASSERT_THROW( PfParam(2, argv4), std::invalid_argument);
    // Check for -eN, -l, and so on ...
    
    //char *argv3[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "7", "8", "5", "-T", "-M", "5.0" };
    //Param pars2 = Param(13, argv3);
    //CPPUNIT_ASSERT_NO_THROW( pars2.parse(model) ); 
    //CPPUNIT_ASSERT_EQUAL( (size_t)3, model.population_number() );
    //CPPUNIT_ASSERT_EQUAL( model.sample_population(4), (size_t)0 );
    //CPPUNIT_ASSERT_EQUAL( model.sample_population(10), (size_t)1 );
    //CPPUNIT_ASSERT_EQUAL( model.sample_population(17), (size_t)2 );
    //CPPUNIT_ASSERT_EQUAL( model.sample_time(4), (double)0.0 );
    //CPPUNIT_ASSERT_EQUAL( model.sample_time(17), (double)0.0 );
    //CPPUNIT_ASSERT_EQUAL( true, pars2.tree_bool );

    //char *argv32[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "7", "8", "5", "5.0", "-T"};
    //CPPUNIT_ASSERT_NO_THROW( Param(12, argv32).parse(model) ); 
    //CPPUNIT_ASSERT_EQUAL( (size_t)3, model.population_number() );
    //CPPUNIT_ASSERT_EQUAL( model.sample_population(4), (size_t)0 );
    //CPPUNIT_ASSERT_EQUAL( model.sample_population(10), (size_t)1 );
    //CPPUNIT_ASSERT_EQUAL( model.sample_population(17), (size_t)2 );
    //CPPUNIT_ASSERT_EQUAL( model.sample_time(4), (double)0.0 );
    //CPPUNIT_ASSERT_EQUAL( model.sample_time(17), (double)0.0 );
    //CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(1, 0)) );
    //CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    //CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 2)) );
    //CPPUNIT_ASSERT_EQUAL( true, pars2.tree_bool );

    //char *argv33[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "7", "8", "5", "5.0"};
    //CPPUNIT_ASSERT_NO_THROW( Param(11, argv33).parse(model) ); 
    //CPPUNIT_ASSERT_EQUAL( (size_t)3, model.population_number() );
    //CPPUNIT_ASSERT_EQUAL( model.sample_population(4), (size_t)0 );
    //CPPUNIT_ASSERT_EQUAL( model.sample_population(10), (size_t)1 );
    //CPPUNIT_ASSERT_EQUAL( model.sample_population(17), (size_t)2 );
    //CPPUNIT_ASSERT_EQUAL( model.sample_time(4), (double)0.0 );
    //CPPUNIT_ASSERT_EQUAL( model.sample_time(17), (double)0.0 );
    //CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(1, 0)) );
    //CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    //CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 2)) );

    //char *argv4[] = { "scrm", "23", "10", "-t", "3.74", "-I", "3", "7", "8", "5", "-eI", "12.3", "2", "0", "1" , "-M", "5.0" };
    //CPPUNIT_ASSERT_NO_THROW( Param(17, argv4).parse(model) ); 
    //CPPUNIT_ASSERT_EQUAL( model.sample_population(20), (size_t)0 );
    //CPPUNIT_ASSERT_EQUAL( model.sample_population(22), (size_t)2 );
    //CPPUNIT_ASSERT_EQUAL( model.sample_time(20), (double)12.3 );
    //CPPUNIT_ASSERT_EQUAL( model.sample_time(22), (double)12.3 );
   
    //// -N & -eN 
    //char *argv5[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "7", "8", "5", 
                      //"-N", "0.3", "-eN", "8.2", "0.75", "-G", "1.5", "-M", "5.0"};
    //CPPUNIT_ASSERT_NO_THROW( Param(19, argv5).parse(model) ); 
    //model.resetTime();
    //CPPUNIT_ASSERT_EQUAL( 0.0, model.getCurrentTime() );
    //CPPUNIT_ASSERT( areSame(0.3*model.default_pop_size, model.population_size(0)) );
    //CPPUNIT_ASSERT( areSame(0.3*model.default_pop_size, model.population_size(1)) );
    //CPPUNIT_ASSERT( areSame(0.3*model.default_pop_size, model.population_size(2)) );
    //model.increaseTime();
    //CPPUNIT_ASSERT( areSame(8.2 * 4 * model.default_pop_size, model.getCurrentTime()) );
    //CPPUNIT_ASSERT( areSame(0.75*model.default_pop_size, model.population_size(0)) );
    //CPPUNIT_ASSERT( areSame(0.75*model.default_pop_size, model.population_size(1)) );
    //CPPUNIT_ASSERT( areSame(0.75*model.default_pop_size, model.population_size(2)) );
    //CPPUNIT_ASSERT_EQUAL( 0.0, model.growth_rate(0) );
    //CPPUNIT_ASSERT_EQUAL( 0.0, model.growth_rate(1) );
    //CPPUNIT_ASSERT_EQUAL( 0.0, model.growth_rate(2) );

    //// -n & -en 
    //char *argv6[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "7", "8", "5", "-G", "1.5", 
                      //"-n", "2", "0.3", "-eN", "1.1", "0.75", "-en", "2", "3", "0.1", "-eG", "1.5", "2", "-M", "5.0" };
    //CPPUNIT_ASSERT_NO_THROW( Param(27, argv6).parse(model) ); 
    //model.finalize();
    ////std::cout << model << std::endl;
    //model.resetTime();
    //CPPUNIT_ASSERT_EQUAL( 0.0, model.getCurrentTime() );
    //CPPUNIT_ASSERT( areSame(model.default_pop_size, model.population_size(0)) );
    //CPPUNIT_ASSERT( areSame(0.3*model.default_pop_size, model.population_size(1)) );
    //CPPUNIT_ASSERT( areSame(model.default_pop_size, (size_t)model.population_size(2)) );
    //CPPUNIT_ASSERT_EQUAL( 1.5, model.growth_rate(0) );
    //CPPUNIT_ASSERT_EQUAL( 0.0, model.growth_rate(1) );
    //CPPUNIT_ASSERT_EQUAL( 1.5, model.growth_rate(2) );
    //model.increaseTime();
    //CPPUNIT_ASSERT_EQUAL( 1.1 * 4 * model.default_pop_size, model.getCurrentTime() );
    //model.increaseTime();
    //CPPUNIT_ASSERT_EQUAL( 1.5 * 4 * model.default_pop_size, model.getCurrentTime() );
    //model.increaseTime();
    //CPPUNIT_ASSERT_EQUAL( 2.0 * 4 * model.default_pop_size, model.getCurrentTime() );
    //CPPUNIT_ASSERT( 0.75*model.default_pop_size > model.population_size(0) );
    //CPPUNIT_ASSERT( 0.75*model.default_pop_size > model.population_size(1) );
    //CPPUNIT_ASSERT_EQUAL( (size_t)(0.10*model.default_pop_size), (size_t)model.population_size(2) );
    //CPPUNIT_ASSERT_EQUAL( 2.0, model.growth_rate(0) );
    //CPPUNIT_ASSERT_EQUAL( 2.0, model.growth_rate(1) );
    //CPPUNIT_ASSERT_EQUAL( 0.0, model.growth_rate(2) );
  }

    void test_initialize_model(){
        int default_nsam = 2;
        double default_mut_rate = 0.001;
        double default_recomb_rate = 0.00001;
        double default_loci_length = 50000;  
    
        char *argv1[] = {"smcsmc"};
        PfParam pfARG_para1(1, argv1);
        vcf_file =  new VariantReader(pfARG_para1.input_variantFileName, VCF, pfARG_para1.buff_length);
        
        Param * scrm_para = new Param(1, argv1, false);    
        CPPUNIT_ASSERT_NO_THROW( scrm_para->parse(model) );
        
        Model * model1 = new Model();        
        //initialize_model(model1, scrm_para, vcf_file, default_nsam, default_mut_rate, default_recomb_rate, default_loci_length);
        CPPUNIT_ASSERT_EQUAL( (size_t)default_nsam, model1->sample_size());
        CPPUNIT_ASSERT_EQUAL( (double)default_mut_rate, model1->mutation_rate());
        CPPUNIT_ASSERT_EQUAL( (double)default_recomb_rate, model1->recombination_rate());
        CPPUNIT_ASSERT_EQUAL( (size_t)default_loci_length, model1->loci_length());
        
        delete vcf_file;
        delete model1;
        delete scrm_para;
        
        char *argv2[] = { "smcsmc", "-vcf", "test6sample.vcf"};
        PfParam pfARG_para2(3, argv2);
        vcf_file =  new VariantReader(pfARG_para2.input_variantFileName, VCF, pfARG_para2.buff_length);
        CPPUNIT_ASSERT_EQUAL( (size_t)3, vcf_file->nsam());
        
        Param * scrm_para2 = new Param(3, argv2, false);    
        CPPUNIT_ASSERT_NO_THROW( scrm_para2->parse(model) );
        
        Model * model2 = new Model();        
        //initialize_model(model2, scrm_para2, vcf_file, default_nsam, default_mut_rate, default_recomb_rate, default_loci_length);
        
        CPPUNIT_ASSERT_EQUAL( (size_t)6, model2->sample_size());
        CPPUNIT_ASSERT_EQUAL( (double)default_mut_rate, model2->mutation_rate());
        CPPUNIT_ASSERT_EQUAL( (double)default_recomb_rate, model2->recombination_rate());
        CPPUNIT_ASSERT_EQUAL( (size_t)default_loci_length, model2->loci_length());
        
        delete vcf_file;
        delete model2;
        delete scrm_para2;
        
    }


  //void testParseMigrationOptions() {
    //// -ma
    //char *argv[] = { "scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
                      //"-ma", "x", "5", "7", "x" };
    //CPPUNIT_ASSERT_NO_THROW( Param(14, argv).parse(model); );
    //model.resetTime();
    //CPPUNIT_ASSERT( areSame(5.0/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    //CPPUNIT_ASSERT( areSame(7.0/(4*model.default_pop_size), model.migration_rate(1, 0)) );

    //// -ema
    //char *argv2[] = { "scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
                      //"-ema", "1.6", "x", "5", "7", "x" };
    //CPPUNIT_ASSERT_NO_THROW( Param(15, argv2).parse(model); );
    //model.resetTime();
    //model.increaseTime();
    //CPPUNIT_ASSERT_EQUAL( 1.6 * 4 * model.default_pop_size, model.getCurrentTime() );
    //CPPUNIT_ASSERT( areSame(5.0/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    //CPPUNIT_ASSERT( areSame(7.0/(4*model.default_pop_size), model.migration_rate(1, 0)) );

    //// -M
    //char *argv3[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "10", "10", "0", 
                      //"-M", "5" };
    //CPPUNIT_ASSERT_NO_THROW( Param(12, argv3).parse(model); );
    //model.resetTime();
    //CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(1, 0)) );
    //CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    //CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 2)) );

    //// -eM
    //char *argv4[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "10", "10", "0", 
                      //"-eM", "1.6", "5" };
    //CPPUNIT_ASSERT_NO_THROW( Param(13, argv4).parse(model); );
    //model.resetTime();
    //model.increaseTime();
    //CPPUNIT_ASSERT_EQUAL( 1.6 * 4 * model.default_pop_size,  model.getCurrentTime() );
    //CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(1, 0)) );
    //CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    //CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 2)) );

    //// -es
    //char *argv5[] = { "scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
                      //"-es", "1.6", "2", "1", "0.5", "-es", "1.6", "1", "2", "0.1", "-M", "5.0" };
    //CPPUNIT_ASSERT_NO_THROW( Param(21, argv5).parse(model); );
    //model.resetTime();
    //model.increaseTime();
    //CPPUNIT_ASSERT_EQUAL( 1.6 * 4 * model.default_pop_size, model.getCurrentTime() );
    //CPPUNIT_ASSERT_EQUAL( 0.1, model.single_mig_pop(0, 1) );
    //CPPUNIT_ASSERT_EQUAL( 0.5, model.single_mig_pop(1, 0) );

    //// -ej
    //char *argv6[] = { "scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
                      //"-ej", "1.6", "2", "1", "-M", "1.3" };
    //CPPUNIT_ASSERT_NO_THROW( Param(15, argv6).parse(model); );
    //model.resetTime();
    //model.increaseTime();
    //CPPUNIT_ASSERT_EQUAL( 1.6 * 4 * model.default_pop_size, model.getCurrentTime() );
    //CPPUNIT_ASSERT_EQUAL( 1.0, model.single_mig_pop(1, 0) );
    //CPPUNIT_ASSERT_EQUAL( 0.0, model.single_mig_pop(0, 1) );

    //CPPUNIT_ASSERT_EQUAL( 0.0, model.migration_rate(0, 1) );
    //CPPUNIT_ASSERT_EQUAL( 1.3/(4*model.default_pop_size), model.migration_rate(1, 0) );
  //}

  void testReadInput() {
    //CPPUNIT_ASSERT_EQUAL( (int)1, readInput<int>("1") );
    //CPPUNIT_ASSERT_EQUAL( (size_t)7, readInput<size_t>("7") );
    //CPPUNIT_ASSERT_EQUAL( (double)3.1, readInput<double>("3.1") );
    //CPPUNIT_ASSERT_THROW( readInput<int>("ABC"), boost::bad_lexical_cast );
    //CPPUNIT_ASSERT_THROW( readInput<int>("-I"), boost::bad_lexical_cast );
  }

  //void testParseGrowthOptions() {
    //// -G && -eG
    //char *argv[] = { "scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
      //"-G", "2", "-eG", "1", "3", "-M", "5.0" };
    //CPPUNIT_ASSERT_NO_THROW( Param(16, argv).parse(model); );
    //model.resetTime();
    //CPPUNIT_ASSERT_EQUAL( 2.0, model.growth_rate(0) );
    //CPPUNIT_ASSERT_EQUAL( 2.0, model.growth_rate(1) );
    //model.increaseTime();
    //CPPUNIT_ASSERT_EQUAL( 1.0 * 4 * model.default_pop_size, model.getCurrentTime() );
    //CPPUNIT_ASSERT_EQUAL( 3.0, model.growth_rate(0) );
    //CPPUNIT_ASSERT_EQUAL( 3.0, model.growth_rate(1) );

    // -g && -eg
    //char *argv2[] = { 
      //"scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
      //"-g", "2", "0.1", "-eG", "1", "3", "-eg", "2", "1", "2.4", "-M", "5.0"};
    //CPPUNIT_ASSERT_NO_THROW( Param(21, argv2).parse(model); );
    //model.resetTime();
    //CPPUNIT_ASSERT_EQUAL( model.default_growth_rate, model.growth_rate(0) );
    //CPPUNIT_ASSERT_EQUAL( 0.1, model.growth_rate(1) );
    //model.increaseTime();
    //CPPUNIT_ASSERT_EQUAL( 1.0 * 4 * model.default_pop_size, model.getCurrentTime() );
    //model.increaseTime();
    //CPPUNIT_ASSERT_EQUAL( 2.0 * 4 * model.default_pop_size, model.getCurrentTime() );
    //CPPUNIT_ASSERT( areSame(2.4, model.growth_rate(0)) );
    //CPPUNIT_ASSERT( areSame(3.0, model.growth_rate(1)) );
  //}
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestParam );
