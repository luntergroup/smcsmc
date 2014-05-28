/*
 * pf-ARG is short for particle filters for ancestral recombination graphs. 
 * This is a free software for demographic inference from genome data with particle filters. 
 * 
 * Copyright (C) 2013, 2014 Sha (Joe) Zhu and Gerton Lunter
 * 
 * This file is part of pf-ARG.
 * 
 * pf-ARG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/


#include"vcf.hpp"
using namespace std;

/*! Initialize vcf file, search for the end of the vcf header. 
 *  Extract the first block of data ( "buffer_length" lines ) into buff
 */
Vcf::Vcf(string file_name, int buffer_length){ /*! Initialize by read in the vcf header file */
    this->init(file_name, buffer_length);
    
    ifstream in_file;
    in_file.open(file_name.c_str());
    in_file.seekg (0, in_file.end);
    vcf_length_ = in_file.tellg();
    
    in_file.seekg (0, in_file.beg);
    string line;
    header_end_pos_=0;
    header_end_line=0;
    in_file.seekg (0, in_file.end);
    vcf_file_length = in_file.tellg();
    in_file.seekg (0, in_file.beg);
    if (in_file.good()){
        getline (in_file,line);
        header_end_pos_ += line.size()+1;
        while (line.size()>0 ){   
            //dout << header_end_line<<"  " <<line.size() <<"  " << header_end_pos_<<"  " << line<<endl;
            if (line[0]=='#'){
                if (line[1]=='#'){
                    // read to header buffer
                    // TO COME ...
                } else {
                    check_feilds(line);
                    //break; //end of the header
                }
            }
    
            getline (in_file,line);
            
            if (line.find("PASS")!= std::string::npos){
                break;
            } // This will skip to the first line of the useable data
            
            header_end_pos_ += line.size()+1;
            header_end_line++;
        }    
    
    }
    this->end_pos_ = this->header_end_pos_;
    in_file.close();
    this->read_new_block(); // END by start reading a block
    }
    

void Vcf::reset_VCF_to_data(){
    /*! Reset the data to the first line, end of the header file 
     *  Extract new block of data
     */ 
    this->end_pos_            = this->header_end_pos_;
    this->end_data_           = false;
    this->eof_                = false;
    this->current_line_index_ = 0;    
    this->current_block_line_ = 0;
    this->site_               = 0;
    this->chrom_              = 0;
    this->previous_site_at_   = 0;
    if ( !this->withdata() ){ this->empty_file_line_counter_ = 0;  }
    //this->empty_file_line_counter_ = 0;
    this->read_new_block();
    }
    
    
void Vcf::init(string infile_name, int buffer_length){
    /*! Initialize the Vcf class members
     */ 
    this->current_line_index_        = 0;
    this->empty_file_line_counter_   = 0;
    this->nsam_                      = 0;
    this->nfield_                    = 0;    
    this->current_block_line_        = 0;
    this->site_                      = 0;
    this->chrom_                     = 0;
    this->previous_site_at_          = 0;
    this->vcf_length_                = 0;
    this->even_interval_             = 0.0;
    this->ghost_num_mut            = 0;
    this->eof_                       = false;
    this->end_data_                  = false;
    this->buffer_max_number_of_lines = buffer_length;
    this->file_name_                 = infile_name;    
    this->withdata_ = ( file_name_.size() > 0 ) ? true : false;    
    }


void Vcf::read_new_line(){
    /*! Read Vcf data, extract mutation site and haplotype
     */ 
    current_line_index_++;
    
    alt.clear();
    vec_of_sample_alt.clear();
    vec_of_sample_alt_bool.clear();
    sample_alt.clear();
    phased.clear(); // True if it is phased, which has '|'
    
    if (!withdata()){
        this->empty_file_line_counter_ ++; 
        this->site_ = ( current_line_index_ == 1 ) ? 0 : this->site_ + even_interval_;
        // add counter for empty file, first time calling read_new_line, site() = 0, second time, set site() = 100000, and haplotype  =  false  
        this->end_data_ = (this->empty_file_line_counter_ > (this->ghost_num_mut +1));
        return;    
    }
    
    string line = this->buffer_lines[current_block_line_];
    skip=true;
    size_t feild_start=0;
    size_t feild_end=0;
    int counter=0;
    string tmp_str;
    while(feild_end<line.size()){
        feild_end=min(line.find('\t',feild_start),line.find('\n',feild_start)); 
        tmp_str=line.substr(feild_start,feild_end-feild_start);
        //istringstream tmp_strm(tmp_str.c_str());
        switch(counter){
            case 0: {
                //tmp_strm>>chrom_; 
                chrom_ = strtod( tmp_str.c_str(), NULL);
                }
            break;
            
            case 1: {
                //tmp_strm>>site_; 
                site_ = strtod( tmp_str.c_str(), NULL);
                if (((site_ - previous_site_at_)<2) && (previous_site_at_ > 0) && (pervious_chrom_ == chrom_)){
                    cout << "Skip reads at chrom " << chrom_<<" at position " <<  site_<<", due to it's too close to the previous variant (at " << previous_site_at_ <<")." << endl;
                    previous_site_at_ = site_;
                    current_block_line_++;
                    return; 
                    }
                }
            break;

            case 3: {
                ref=tmp_str; 
                //if (ref.size()>1){
                //cout << "Skip reads at chrom " << chrom_<<" at position " <<  site_<<", due to deletion or replacement" << endl;
                //return;}
                }
            break;
            
            case 4: { 
                size_t alt_start=0;size_t alt_end=0; string alt_str;
                while (alt_end<tmp_str.size()){
                    alt_end=min(tmp_str.find(',',alt_start),tmp_str.size());
                    alt_str=tmp_str.substr(alt_start,alt_end);
                    //if (alt_str.size()>1){
                        //cout << "Skip reads at chrom " << chrom_<<" at position " <<  site_<<", due to insertion" << endl;
                        //return;}
                    alt.push_back(alt_str);
                    alt_start=alt_end+1;
                    }
                }
            break;
            
            case 6: {
                if (tmp_str!="PASS"){
                    cout << "Skip reads at chrom " << chrom_<<" at position " <<  site_<<", due to low qualitiy." << endl;
                    current_block_line_++;
                    return;
                    }
                }
            break;
            }

        if (counter > 8){
            size_t bar_index=tmp_str.find('|',0);        
            size_t slash_index=tmp_str.find('/',0);
            size_t colon_index=tmp_str.find(':',0);
            size_t break_index=min(bar_index, slash_index);
            assert(break_index<colon_index);
            vec_of_sample_alt.push_back(extract_alt_(tmp_str, 0, break_index));            
            vec_of_sample_alt.push_back(extract_alt_(tmp_str, break_index+1, colon_index));            
            size_t alt_index_0 = strtol (tmp_str.substr(0,1).c_str(), NULL, 0);;
            vec_of_sample_alt_bool.push_back( (alt_index_0 == (size_t)0) ? false : true);
            size_t alt_index_2 = strtol (tmp_str.substr(2,1).c_str(), NULL, 0);;
            vec_of_sample_alt_bool.push_back( (alt_index_2 == (size_t)0) ? false : true);
            }
        feild_start=feild_end+1;        
        counter++;
        }
    skip=false;
    // CHECK END of Vcf
    
    current_block_line_++;
    if (current_block_line_ == buffer_lines.size() ){ // END of the buff block
        if (!this->eof_){
            this->empty_block();
            this->read_new_block();
            }
        else{
            this->end_data_=true;    
            }
        }
    }


void Vcf::empty_block() { buffer_lines.clear(); }


void Vcf::read_new_block(){
    
    if (current_line_index_ == 0){
        cout << "Set data to the first entry, read a block of " <<  this->buffer_max_number_of_lines << " entries" <<endl;
        } 
    else {
        cout << "Read new block of " <<  this->buffer_max_number_of_lines << " entries, currently at line " << current_line_index_ << endl;
        }
    
    buffer_lines.clear();
    current_block_line_ = 0;

    if (!this->withdata_){
        this->site_ = 0;
        this->eof_=true;
        return;
        }

    ifstream infile (file_name_.c_str(), std::ifstream::binary);
    infile.seekg(this->end_pos_, ios::beg);
    int count=0;
    string tmp_str;
    
    getline(infile, tmp_str);
    
    while (infile.good() && tmp_str.size()>0 && count < buffer_max_number_of_lines){
        buffer_lines.push_back(tmp_str);
        this->end_pos_ += tmp_str.size()+1;
        count++;
        getline(infile, tmp_str);
        }
    
    this->eof_ = ( vcf_length_ <= this->end_pos_ );
    
    infile.close();
    }


string Vcf::extract_alt_(string tmp_str, size_t start, size_t end){
    /*! Extract haplotype
     */ 
    size_t alt_index = strtol (tmp_str.substr(start,end-start).c_str(), NULL, 0);
    string alt_dummy = ( alt_index==0 ) ? ref : alt[alt_index-1];
    return alt_dummy;
    }
    

void Vcf::check_feilds(string line){
    size_t feild_start=0;
    size_t feild_end=0;
    int counter=0;
    string tmp_str;
    while( feild_end < line.size() ) {
        feild_end = min( line.find('\t',feild_start), line.find('\n',feild_start) ); 
        tmp_str = line.substr( feild_start, feild_end-feild_start );
        switch (counter){
            case 0: if ( tmp_str != "#CHROM" ){ throw std::invalid_argument( "First Header entry should be #CHROM: "   + tmp_str ); } break;
            case 1: if ( tmp_str != "POS"    ){ throw std::invalid_argument( "Second Header entry should be POS: "     + tmp_str ); } break;
            case 2: if ( tmp_str != "ID"     ){ throw std::invalid_argument( "Third Header entry should be ID: "       + tmp_str ); } break;
            case 3: if ( tmp_str != "REF"    ){ throw std::invalid_argument( "Fourth Header entry should be REF: "     + tmp_str ); } break; 
            case 4: if ( tmp_str != "ALT"    ){ throw std::invalid_argument( "Fifth Header entry should be ALT: "      + tmp_str ); } break;
            case 5: if ( tmp_str != "QUAL"   ){ throw std::invalid_argument( "Sixth Header entry should be QUAL: "     + tmp_str ); } break;
            case 6: if ( tmp_str != "FILTER" ){ throw std::invalid_argument( "Seventh Header entry should be FILTER: " + tmp_str ); } break;
            case 7: if ( tmp_str != "INFO"   ){ throw std::invalid_argument( "Eighth Header entry should be INFO: "    + tmp_str ); } break;
            case 8: if ( tmp_str != "FORMAT" ){ throw std::invalid_argument( "Ninth Header entry should be FORMAT: "   + tmp_str ); } break;
            }
        
        if (counter > 8){  sample_names.push_back(tmp_str); }
        
        feild_start = feild_end+1;
        counter++;
        } // End of while loop: feild_end < line.size()
    
    this->nfield_ = counter;
    set_nsam( (int)sample_names.size() );
    assert( print_sample_name() );   
    
    }



