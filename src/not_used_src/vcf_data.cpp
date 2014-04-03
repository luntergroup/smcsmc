//vcf_data.cpp 
#include"vcf_data.hpp"
using namespace std;
using namespace vcf_data;


vcf_header::vcf_header(string file_name, int usr_buffer_number_of_lines){
	buffer_max_number_of_lines=usr_buffer_number_of_lines;
	ifstream in_file;
	in_file.open(file_name.c_str());
	string line;
	header_end_pos=0;
	header_end_line=0;
	in_file.seekg (0, in_file.end);
    vcf_file_length = in_file.tellg();
    in_file.seekg (0, in_file.beg);
	if (in_file.good()){
		getline (in_file,line);
		header_end_pos=header_end_pos+line.size()+1;
		//while (line.size()>0){   
		while (line.size()>0 ){   
		dout<<header_end_line<<"  " <<line.size() <<"  " << header_end_pos<<"  " << line<<endl;
			if (line[0]=='#'){
				if (line[1]=='#'){
					// read to header buffer
					// TO COME ...
				}
				else{
					check_feilds(line);
					//break; //end of the header
				}
			}
			
			getline (in_file,line);
			if (line.find("PASS")!= std::string::npos){break;}
			header_end_pos=header_end_pos+line.size()+1;
			header_end_line++;
		}	
		
	}
	
	//in_file.seekg (0, in_file.end);
    //max_buff_length = in_file.tellg() ;
    //max_buff_length= max_buff_length - int(header_end_line);
    //cout << "max_buff_length =" << max_buff_length<<endl;
	in_file.close();
}

vcf_buffer::vcf_buffer(string file_name, size_t previous_end_pos, int usr_buffer_number_of_lines){
	if (file_name.size()==0){
		
		string tmp_str="#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001";
		buffer_lines.push_back(tmp_str);
		tmp_str="1	1	rs6040355	A	T	67	PASS	NS=2;	GT	0|0";
		buffer_lines.push_back(tmp_str);
		tmp_str="1	50000	rs6040355	A	T	67	PASS	NS=2;	GT	0|0";
		buffer_lines.push_back(tmp_str);
		return;
	}
	//cout << "Read Buffer ..." << endl;
	buffer_max_number_of_lines=usr_buffer_number_of_lines;
	//cout << "shoudl read " << buffer_max_number_of_lines<<" lines" << endl;
	end_pos=previous_end_pos;
	ifstream infile (file_name.c_str(), std::ifstream::binary);
	infile.seekg(end_pos,ios::beg);
	int count=0;
	string tmp_str;
	getline(infile, tmp_str);

	while (infile.good() && tmp_str.size()>0 && count < buffer_max_number_of_lines){
	
		cout<<tmp_str<<endl;		
		buffer_lines.push_back(tmp_str);
		end_pos=end_pos+tmp_str.size()+1;
		count++;
		getline(infile, tmp_str);
		//if (tmp_str.find("PASS")!= std::string::npos){
			//break;
		//}
	}
	//cout<<count << " lines" << endl;
    infile.close();
}

void vcf_header::check_feilds(string line){
	size_t feild_start=0;
	size_t feild_end=0;
	int counter=0;
	string tmp_str;
	while(feild_end<line.size()){
		feild_end=min(line.find('\t',feild_start),line.find('\n',feild_start)); 
		tmp_str=line.substr(feild_start,feild_end-feild_start);
		switch (counter){
			case 0: if (tmp_str != "#CHROM"){ throw std::invalid_argument("First Header entry should be #CHROM: " + tmp_str);} break;
			case 1: if (tmp_str != "POS"){ throw std::invalid_argument("Second Header entry should be POS: " + tmp_str);} break;
			case 2: if (tmp_str != "ID"){ throw std::invalid_argument("Third Header entry should be ID: " + tmp_str);} break;
			case 3: if (tmp_str != "REF"){ throw std::invalid_argument("Fourth Header entry should be REF: " + tmp_str);}break; 
			case 4: if (tmp_str != "ALT"){throw std::invalid_argument("Fifth Header entry should be ALT: " + tmp_str);}break;
			case 5: if (tmp_str != "QUAL"){throw std::invalid_argument("Sixth Header entry should be QUAL: " + tmp_str);}break;
			case 6: if (tmp_str != "FILTER"){ throw std::invalid_argument("Seventh Header entry should be FILTER: " + tmp_str);}break;
			case 7: if (tmp_str != "INFO"){ throw std::invalid_argument("Eighth Header entry should be INFO: " + tmp_str);}break;
			case 8: if (tmp_str != "FORMAT"){ throw std::invalid_argument("Ninth Header entry should be FORMAT: " + tmp_str);}break;		
		}
			//cout<<counter<<"  " << tmp_str<<"  " << feild_start << " " <<feild_end <<endl;		
		if (counter > 8){
			sample_names.push_back(tmp_str);
		}
		feild_start=feild_end+1;		
		counter++;
	}
	nfield_=counter;
	set_nsam(int(sample_names.size()));
	dout << "Sample names:" << endl;for (int i=0;i<nsam();i++){	dout<<sample_names[i]<<" ";	}dout<<endl;
}


vcf_line::vcf_line(string line,int nsam,int pervious_chrom,double previous_site_at, int first_sample_field, int poistion_field){
	cout<<line<<endl;
	skip=true;
	size_t feild_start=0;
	size_t feild_end=0;
	int counter=0;
	string tmp_str;
	while(feild_end<line.size()){
		feild_end=min(line.find('\t',feild_start),line.find('\n',feild_start)); 
		tmp_str=line.substr(feild_start,feild_end-feild_start);
		istringstream tmp_strm(tmp_str.c_str());
		switch(counter){
			case 0: tmp_strm>>chrom_; break;
			case 1: {tmp_strm>>site_; 
				//if (((site_ - previous_site_at)<150) && (previous_site_at > 0) && (pervious_chrom==chrom_)){
				if (((site_ - previous_site_at)<2) && (previous_site_at > 0) && (pervious_chrom==chrom_)){
					cout << "Skip reads at chrom " << chrom_<<" at position " <<  site_<<", due to it's too close to the previous variant (at " << previous_site_at <<")." << endl;
					return; // in nsam = 0, then do not skip ...
					}
			}break;
			case 3: {ref=tmp_str; 
				//if (ref.size()>1){
					//cout << "Skip reads at chrom " << chrom_<<" at position " <<  site_<<", due to deletion or replacement" << endl;
					//return;}
				} break;
			case 4: { size_t alt_start=0;size_t alt_end=0; string alt_str;
				while (alt_end<tmp_str.size()){
					alt_end=min(tmp_str.find(',',alt_start),tmp_str.size());
					alt_str=tmp_str.substr(alt_start,alt_end);
					//if (alt_str.size()>1){
						//cout << "Skip reads at chrom " << chrom_<<" at position " <<  site_<<", due to insertion" << endl;
						//return;}
					alt.push_back(alt_str);
					alt_start=alt_end+1;
				}
				}break;
			case 6: {
				if (tmp_str!="PASS"){
					cout << "Skip reads at chrom " << chrom_<<" at position " <<  site_<<", due to low qualitiy." << endl;
					return;}
				}break;
		}

		if (counter > 8){
			size_t bar_index=tmp_str.find('|',0);		
			size_t slash_index=tmp_str.find('/',0);
			size_t colon_index=tmp_str.find(':',0);
			size_t break_index=min(bar_index, slash_index);
			assert(break_index<colon_index);
			vec_of_sample_alt.push_back(extract_alt_(tmp_str, 0, break_index));
			
			vec_of_sample_alt.push_back(extract_alt_(tmp_str, break_index+1, colon_index));
			
			istringstream alt_index_0_strm(tmp_str.substr(0,1));
			size_t alt_index_0;
			alt_index_0_strm>>alt_index_0;
			if (alt_index_0==0){vec_of_sample_alt_bool.push_back(false);}else{vec_of_sample_alt_bool.push_back(true);}
			
			istringstream alt_index_2_strm(tmp_str.substr(2,1));
			size_t alt_index_2;
			alt_index_2_strm>>alt_index_2;
			if (alt_index_2==0){vec_of_sample_alt_bool.push_back(false);}else{vec_of_sample_alt_bool.push_back(true);}

			
			////consider to remove the following ..
			//vector <string> sample_alt_dummy;
			//sample_alt_dummy.push_back(extract_alt_(tmp_str, 0, break_index));
			//sample_alt_dummy.push_back(extract_alt_(tmp_str, break_index+1, colon_index));
			//sample_alt.push_back(sample_alt_dummy);
			////consider to remove ...
			
			
		}
		feild_start=feild_end+1;		
		counter++;
	}
	skip=false;
}


//bool vcf_line::extract_alt_bool_(string tmp_str, size_t start, size_t end){
	//string alt_index_str=tmp_str.substr(start,end-start);
	//istringstream alt_index_strm(alt_index_str.c_str());
	//size_t alt_index;
	//alt_index_strm>>alt_index;

	//bool alt_dummy;
	//if (alt_index==0){
		//alt_dummy=ref;
	//}
	//else{
		//alt_dummy=alt[alt_index-1];
	//}
	//return alt_dummy;
//}


string vcf_line::extract_alt_(string tmp_str, size_t start, size_t end){
	string alt_index_str=tmp_str.substr(start,end-start);
	//cout<<alt_index_str<<endl;
	//cout << "hea" << endl;
	istringstream alt_index_strm(alt_index_str.c_str());
	size_t alt_index;
	alt_index_strm>>alt_index;

	string alt_dummy;
	if (alt_index==0){
		alt_dummy=ref;
	}
	else{
		alt_dummy=alt[alt_index-1];
	}
	//cout << "extract_alt_ finished" << endl;
	return alt_dummy;
}

void vcf_line::print_vcf_line(vector<string> sample_names){
	dout << "CHROM = " <<chrom_<<endl;
	dout << "POS = " << site_ <<endl;
	for (size_t i=0;i<sample_names.size();i++){
		dout<<sample_names[i]<<": ";
		for (size_t j=0;j<sample_alt[i].size();j++){
			dout<<sample_alt[i][j]<<" ";
		}
		cout<<endl;
	}
	dout<<endl;
}
