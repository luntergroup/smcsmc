#include"vcf.hpp"
#include"pfparam.hpp"
#include"usage.hpp"


using namespace std;



int main(int argc, char *argv[]){
    /*!
     * Settings 
     */
    //string defaultvcf("testdata.vcf");
    string defaultvcf("../pfpsmc/R/sim-YHdata2.vcf");
    int default_buff_length = 100;
    
    if (argc > 1){
        string argv_i(argv[1]);
        defaultvcf = argv_i;
    }
    
    
    //INITIALIZE USAGE
    int who = RUSAGE_SELF;
    struct rusage usage;
    struct rusage *p=&usage;


    /*! 
     * Start testing, expect to see no grow in memory when load data
     */
      
    pfARG::param pfARG_para(argc, argv);
    pfARG_para.vcf_NAME = defaultvcf;
    pfARG_para.buff_length = default_buff_length;
    
    
    Vcf * VCFfile =  new Vcf(pfARG_para.vcf_NAME, pfARG_para.buff_length);
    
    cout << pfARG_para.vcf_NAME << endl;
    
    VCFfile->read_new_line();
    int counter = 1;
    cout << VCFfile->nsam() << endl;
    do{ 
        VCFfile->read_new_line(); // Read new line from the vcf file      
        counter++;
        printCurrentUsage(who, p, 0, 0); //check for the usage        

    }while(!VCFfile->end_data());
    
    delete VCFfile;
    cout << "max buff size = " << default_buff_length << endl;
    return 0;    
}
