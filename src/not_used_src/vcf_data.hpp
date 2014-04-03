//vcf_data.hpp
#include"global.hpp"
#include<string>
//#include <iostream>
#include <sstream>
using namespace std;

namespace vcf_data{
	class vcf_header{
		public:
			vcf_header(string file_name, int usr_buffer_max_of_lines=3);
			~vcf_header(){};
			int nsam() const {return this->nsam_;};
			void set_nsam(int nsam) {this->nsam_ = nsam;};
			
			int nfield() const{return this->nfield_;};	
			size_t vcf_file_length;
			vector <string> sample_names;
			//void scan_more();
			//char ** buffer;
			
			size_t header_end_line;
			size_t header_end_pos;
			//int buffer_start_line;
			int buffer_max_number_of_lines;
			//int max_buff_length;
		private:
			int nsam_;
			//bool eof;
			//void read_header(string line);
			
			void check_feilds(string line);
			int nfield_;
		
	};
	
	
	class vcf_buffer{
		public:
		vcf_buffer(string file_name, size_t previous_end_pos, int usr_buffer_max_of_lines=3);
		
		vcf_buffer(){};
		~vcf_buffer(){};
		size_t end_pos;

		int buffer_max_number_of_lines;
		vector <string> buffer_lines;
		private:
		
		
	};
	
	
	class vcf_line{
		public:
			vcf_line(string line,int nsam, int pervious_chrom, double previous_site_at=0,int first_sample_field=9, int poistion_field=1);
			~vcf_line(){};
			
			double site() const{return this->site_;};
			int chrom() const{return this ->chrom_;};
			string ref;
			vector<string> alt;
			vector <string> vec_of_sample_alt;
			vector <bool> vec_of_sample_alt_bool;
			vector < vector<string> > sample_alt;
			vector<bool> phased; // True if it is phased, which has '|'
			bool skip;
			
			void print_vcf_line(vector<string> sample_names);
			
		private:
			string extract_alt_(string tmp_str, size_t start, size_t end);
			//bool extract_alt_bool_(string tmp_str, size_t break_index);
			double site_;
			int chrom_;
	};
	
}
