/*
 * smcsmc is short for particle filters for ancestral recombination graphs. 
 * This is a free software for demographic inference from genome data with particle filters. 
 * 
 * Copyright (C) 2013, 2014 Sha (Joe) Zhu and Gerton Lunter
 * 
 * This file is part of smcsmc.
 * 
 * smcsmc is free software: you can redistribute it and/or modify
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


#include"variantReader.hpp"
using namespace std;

/*! Initialize vcf file, search for the end of the vcf header. 
 *  Extract the first block of data ( "buffer_length" lines ) into buff
 */
VariantReader::VariantReader(string file_name, INPUT_FILETYPE FileType_in, int buffer_length): buffer_max_number_of_lines(buffer_length) { 
    //try{
    this->file_name_ = file_name;
    this->FileType = this->file_name_.size() > 0 ? FileType_in : EMPTY; 
    /*! Initialize by read in the vcf header file */
    this->init();

    if ( this->FileType != EMPTY ){
        this->skip_Header();
    }
    this->read_new_block(); // END by start reading a block
    //}
    //catch (const string &e) {
        ////std::cerr << "Error: " << e.what() << std::endl;
        //std::cerr << "Error: " << e << std::endl;
    //}
}

void VariantReader::skip_Header(){
    ifstream in_file;
    in_file.open( this->file_name_.c_str());
    in_file.seekg (0, in_file.end);
    this->file_length_ = in_file.tellg();
    
    in_file.seekg (0, in_file.beg);
    header_end_pos_=0;
    header_end_line=0;
    //in_file.seekg (0, in_file.end);
    //in_file.seekg (0, in_file.beg);
    if ( in_file.good() ){
        getline ( in_file, this->tmp_line ); 
        header_end_pos_ += this->tmp_line.size() + 1;
        //dout << this->tmp_line <<endl;
        while ( this->tmp_line.size()>0 ){   
            if ( this->tmp_line[0]=='#' ){
                if ( this->tmp_line[1]=='#' ){
                    // read to header buffer
                    // TO COME ...
                } 
                else {
                    this->check_feilds();
                    //break; //end of the header
                }
            }
    
            getline ( in_file, this->tmp_line );
            
            // This requires more thinking, are we checking the data quality?
            if ( this->tmp_line.find("PASS")!= std::string::npos || this->tmp_line.find("REFCALL")!= std::string::npos ){ break; } // This will skip to the first line of the useable data
            
            header_end_pos_ += this->tmp_line.size() + 1;
            header_end_line++;
        }    
    
    }
    else {
        throw std::invalid_argument ( std::string("Invaild file: ") + this->file_name_ );
        }
    this->end_pos_ = this->header_end_pos_;
    in_file.close();
    }
    
void VariantReader::init(){
    /*! Initialize other VariantReader class members
     */
    this->reset_chrom_site(); // reinitilize VariantPosition   
    this->filter_window_             = 1;
    this->current_line_index_        = 0;
    this->empty_file_line_counter_   = 0;
    this->nsam_                      = 0;
    //this->nfield_                    = 0;    
    //this->current_block_line_        = 0;
    this->missing_data_threshold_ = INT_MAX;
    this->file_length_                = 0;
    this->even_interval_             = 0.0;
    this->ghost_num_mut            = 0;
    this->eof_                       = false;
    this->end_data_                  = false;
    this->empty_block();
    }


void VariantReader::reset_data_to_first_entry(){
    /*! Reset the data to the first line, end of the header file 
     *  Extract new block of data
     */
    this->reset_chrom_site(); // reinitilize VariantPosition     
    this->end_pos_            = this->header_end_pos_;
    this->end_data_           = false;
    this->eof_                = false;
    this->current_line_index_ = 0;    
    if ( this->FileType == EMPTY ){ this->empty_file_line_counter_ = 0; }
    this->read_new_block();
    }


void VariantReader::read_new_line(){
    /*! Read VariantReader data, extract mutation site and haplotype
     */ 
    this->current_line_index_++;
    
    if ( this->FileType == EMPTY ){
        this->set_empty_line_entry();
        return;    
    }
    
    this->initialize_read_newLine(); 
    assert ( !this->skip_tmp_line );

    while ( field_end < this->tmp_line.size() && !this->skip_tmp_line ){
        field_end = min ( this->tmp_line.find('\t',feild_start), this->tmp_line.find('\n', feild_start) );
        this->tmp_str = this->tmp_line.substr( feild_start, field_end - feild_start );
        
        if      ( field_index == 0 ) { this->extract_field_CHROM();   }  /*! Throw exception when it is on different chrom*/        
        else if ( field_index == 1 ) { this->extract_field_POS ();    }  /*! Skip the line if two sites are two close, apply the filtering */
        //else if ( field_index == 2 ) { this->extract_field_ID ();   }        
        else if ( field_index == 3 ) { this->extract_field_REF () ;   }  /*! REF field, Skip the line, if the length of REF string is greater than one, treated as deletion or indel*/
        else if ( field_index == 4 ) { this->extract_field_ALT () ;   }  /*! ALT field, Skip the line if any alter field */
                                                                         /*! what about the case of 1|1, ./.*/
        //else if ( field_index == 5 ) { this->extract_field_QUAL() ; } // ATTENTION, DO NOT CHECK QUALITY AT THE MOMENT
        else if ( field_index == 6 ) { this->extract_field_FILTER();  }  /*! Skip if it is not passed for snp, or it is not refcall */
        else if ( field_index == 7 ) { this->extract_field_INFO();    } // READ INFO FIELD, EXTRACT THE LENGTH OF MISSING BASE FROM RGVCF, OR EXTRACT THE LENGTH OF INVARIANT FROM GVCF
        else if ( field_index > 8  ) { this->extract_field_VARIANT(); }
        
        feild_start = field_end+1;        
        field_index++;
        }
    //
    // CHECK END of VariantReader
    this->current_block_line_++;        
    this->check_and_update_block();
    this->check_and_update_newLine();
    this->finalize_read_new_line();
    }

void VariantReader::set_empty_line_entry(){
    this->current_variant_state = INVARIANT ;
    this->previous_seg_state = MISSING;
    this->current_seg_state  = MISSING;

    this->site_ = ( current_line_index_ == 1 ) ? 0 : this->site_ + even_interval_;
    this->seg_end_site_ = this->site_;

    // add counter for empty file, first time calling read_new_line, site() = 0, second time, set site() = 100000, and haplotype  =  false  
    this->empty_file_line_counter_ ++; 
    this->end_data_ = (this->empty_file_line_counter_ > (this->ghost_num_mut + 1)); // DEBUG 
    }    


void VariantReader::initialize_read_newLine(){
    this->previous_variant_state = this->current_variant_state ;

// \todo !!! need to work on the previous_seg_state
    this->previous_seg_state = this->current_seg_state;
    //this->current_seg_state = SEQ_INVARIANT;
    // Initialize read new line
    alt.clear();
    vec_of_sample_alt.clear();
    int_vec_of_sample_alt.clear();
    sample_alt.clear();
    phased.clear(); // True if it is phased, which has '|'    

    this->feild_start = 0;
    this->field_end   = 0;
    this->field_index    = 0;
    this->skip_tmp_line = false;
    
    this->tmp_line = this->buffer_lines[this->current_block_line_];

    }

/*! Check if it has reached the end of the buffer block, if so, read the next block */ 
void VariantReader::check_and_update_block(){
    if ( (this->current_block_line_ == this->buffer_lines.size() ) && !this->eof_ ){ // END of the buff block
        this->read_new_block();
        }
    this->end_data_ = this->eof_;    
    }

void VariantReader::check_and_update_newLine(){
    if ( this->skip_tmp_line ){
        // revert, if the current line is going to be skipped.
        
        this->current_seg_state = this->previous_seg_state ;
        this->current_variant_state = this->previous_variant_state ;
        if ( !this->end_data_ ){ this->read_new_line(); }  // If this line is bad, and hasn't reached to the end of the data. skip, then read again
        else { this->site_ = this->previous_site_at_; }
        }    
    }

void VariantReader::finalize_read_new_line(){

    // If we have made to here, that means, this line is good, so, we need to check the distance between this line and the next line. 
    if ( (site_ - previous_site_at_) > missing_data_threshold_ ){
        cout << " New data ( "<< site_ <<" ) is too far away from previous (" << previous_site_at_ << "), treat as missing data" << endl;
        //this->set_missding_data ( true );
        this->previous_seg_state = MISSING;
        } 
            
    // If it is not skipped, overwrite the previous site and chrom
    this->previous_site_at_ = this->site_; // If the current line is valid, update the previous site to the current site.
    this->pervious_chrom_ = this->chrom_ ; // If the current line is valid, update the previous chrom_ to the current chrom_. to check if they are on the same chrom        
    }

void VariantReader::empty_block() {
    this->current_block_line_ = 0;
    this->buffer_lines.clear(); 
    }


void VariantReader::read_new_block(){
    
    if ( this->FileType == EMPTY ){
        this->site_ = 0;
        this->eof_  = true;
        return;
        }
        
    this->empty_block();
    
    if (current_line_index_ == 0){ cout << "Set data to the first entry, read a block of " <<  this->buffer_max_number_of_lines << " entries" <<endl; } 
    else { cout << "Read new block of " <<  this->buffer_max_number_of_lines << " entries, currently at line " << current_line_index_ << endl;  }
    
    ifstream infile (file_name_.c_str(), std::ifstream::binary);
    infile.seekg( this->end_pos_, ios::beg );
    int line_counter = 0;
    
    getline( infile, this->tmp_line );

    while (infile.good() && this->tmp_line.size()>0 && (line_counter) < buffer_max_number_of_lines){ 
        this->buffer_lines.push_back( this->tmp_line );
        this->end_pos_ += this->tmp_line.size()+1;
        line_counter++;
        getline(infile, this->tmp_line);
        }
    
    this->eof_ = ( this->file_length_ <= this->end_pos_ );
    infile.close();
    }


void VariantReader::extract_field_CHROM () { 

    this->chrom_ = strtol( tmp_str.c_str(), NULL, 0); 
    if ( pervious_chrom_ != chrom_ && pervious_chrom_!= 0 ) { throw std::invalid_argument ( "Two different chroms" ); }
    }


void VariantReader::extract_field_POS ( ){ 
    
    this->site_ = strtol( tmp_str.c_str(), NULL, 0); 

    if (  this->pervious_chrom_ > 0) {assert ( this->pervious_chrom_ == this->chrom_ );}
    assert ( this->previous_site_at_ >= 0);

    if ( ( this->site_ - this->previous_site_at_ ) < this->filter_window_   ) { // Making 100, as a filtering process, to screen mutations that are too close
        cout << "Skip reads at chrom " << this->chrom_<<" at position " <<  this->site_ << ", due to it's too close to the previous variant (at " << this->previous_site_at_ <<")." << endl;
        this->skip_tmp_line = true;
        }
    }

    
void VariantReader::extract_field_ID ( ){}


void VariantReader::extract_field_REF ( ){
    /*! Check the REF field, if there REF length is greater than 1, it is indel.
     *  if it is gvcf file, and REF is "." or "N", it is INVARIANT, 
     */ 
    
    this->ref = tmp_str; 
    
    if ( this->ref.size() > 1 ){ // BAD LINE, SKIP EXACTING INFORMATION FOR FIELD
        cout << "Skip reads at chrom " << this->chrom_<<" at position " <<  this->site_<<", due to deletion or replacement" << endl;
        //this->current_variant_state = OTHER_VARIANT; // IGNORE OTHER_VARIANT FOR NOW, JUST SKIP
        this->skip_tmp_line = true;
        return ;
    } 
    
    if ( this->ref == "." || this->ref == "N" ) { 
        if ( this->FileType == GVCF  ) { 
             this->current_variant_state = INVARIANT; 
             this->current_seg_state = SEQ_INVARIANT;
             }
        else if ( this->FileType == RGVCF ) { 
             this->current_variant_state = INVARIANT;
             this->current_seg_state = MISSING; 
             }
        }
    else { // assume besides ".", "N", it can only be "A","T", "C", "G"
        this->current_variant_state = SNP;
        this->current_seg_state = ZERO_SEG;            
        }
}


void VariantReader::extract_field_ALT ( ){
    if ( this->tmp_str == "." || this->tmp_str == "N" ) { 
        if ( this->FileType == GVCF  ) { 
             this->current_variant_state = INVARIANT; 
             this->current_seg_state = SEQ_INVARIANT;
             }
        if ( this->FileType == RGVCF ) { 
             this->current_variant_state = INVARIANT;
             this->current_seg_state = MISSING; 
             }
        return;
    }    
    size_t alt_start = 0; 
    size_t alt_end = 0; 
    string alt_tmp_str;
    bool alt_not_valid = false; // Depends on the number of samples, need to check each ALT string for each sample
    while ( alt_end < this->tmp_str.size()){
        alt_end = min( this->tmp_str.find(',', alt_start), this->tmp_str.size() );
        alt_tmp_str = this->tmp_str.substr(alt_start, alt_end);
        if ( alt_tmp_str.size() > 1 ){ // BAD LINE, SKIP EXACTING INFORMATION FOR FIELD
            cout << "Skip reads at chrom " << this->chrom_<<" at position " <<  this->site_ <<", due to insertion" << endl;
            this->current_variant_state = OTHER_VARIANT;
            alt_not_valid = true;
            break;
            }
        alt.push_back(alt_tmp_str);
        alt_start=alt_end+1;
        }
    if ( alt_not_valid ){
        this->skip_tmp_line = true;
        }           
    }
    
void VariantReader::extract_field_QUAL ( ){}

void VariantReader::extract_field_FILTER ( ){ 
    this->skip_tmp_line = false;
    if ( this->tmp_str.find( "PASS") == std::string::npos && this->tmp_str.find( "REFCALL") == std::string::npos ){ // BAD LINE, SKIP EXACTING INFORMATION FOR FIELD
        cout << "Skip reads at chrom " << chrom_<<" at position " <<  site_<<", due to low qualitiy." << endl;
        //cout << "skipping: "<< line << endl; // DEBUG
        this->skip_tmp_line = true;
        } 
    }
    
void VariantReader::extract_field_INFO ( ){

    if ( this->current_variant_state == SNP ){ 
        assert( this->current_seg_state == ZERO_SEG );
        this->seg_end_site_ = this->site_; 
        }
    else {
        assert( this->tmp_str.find( "END=", 0 ) != std::string::npos );
        this->seg_end_site_ = strtol( tmp_str.substr( (size_t)4 ).c_str(), NULL, 0);        
        
        assert ( (this->FileType == GVCF) || ( this->FileType == RGVCF ) );
        this->current_seg_state = this->FileType == GVCF ? SEQ_INVARIANT : MISSING ;
        //if      ( this->FileType == VCF )  this->current_seg_state = SEQ_INVARIANT;
        //else if ( this->FileType == GVCF )  this->current_seg_state = SEQ_INVARIANT;
        //else    ( this->FileType == RGVCF )  this->current_seg_state = MISSING;
        
        //if ( this->FileType == GVCF ) { this->current_seg_state = SEQ_INVARIANT; }
        //else if ( this->FileType == RGVCF ) { this->current_seg_state = MISSING; }
        }
    }


void VariantReader::extract_field_FORMAT ( ){ }

void VariantReader::extract_field_VARIANT ( ){
    //assert ( this->skip_tmp_line = false );
    if ( this->current_variant_state == INVARIANT ){
        this->int_vec_of_sample_alt.push_back( 0 );
        this->int_vec_of_sample_alt.push_back( 0 );
        return;
    }
    
    size_t bar_index   = this->tmp_str.find('|',0); 
    size_t slash_index = this->tmp_str.find('/',0);
    size_t colon_index = this->tmp_str.find(':',0);
    size_t break_index = min(bar_index, slash_index);
    assert( break_index < colon_index );
    this->vec_of_sample_alt.push_back( extract_field_ALT_str( 0, break_index));            
    this->vec_of_sample_alt.push_back( extract_field_ALT_str( break_index+1, colon_index));
    
    size_t alt_index_0 = strtol (tmp_str.substr(0,1).c_str(), NULL, 0);;    
    this->int_vec_of_sample_alt.push_back( (alt_index_0 == (size_t)0) ? 0 : 1);
    
    size_t alt_index_2 = strtol (tmp_str.substr(2,1).c_str(), NULL, 0);;
    this->int_vec_of_sample_alt.push_back( (alt_index_2 == (size_t)0) ? 0 : 1);    
    }


string VariantReader::extract_field_ALT_str( size_t start, size_t end ){
    /*! Extract haplotype */ 
    size_t alt_index = strtol ( this->tmp_str.substr(start,end-start).c_str(), NULL, 0);
    return ( alt_index==0 ) ? ref : alt[alt_index-1];
}
    

void VariantReader::check_feilds(){
    size_t feild_start = 0;
    size_t field_end = 0;
    int field_index = 0;
    //string tmp_str;
    while( field_end < this->tmp_line.size() ) {
        field_end = min( this->tmp_line.find('\t',feild_start), this->tmp_line.find('\n',feild_start) ); 
        this->tmp_str = this->tmp_line.substr( feild_start, field_end-feild_start );
        switch ( field_index ){
            case 0: if ( this->tmp_str != "#CHROM" ){ throw std::invalid_argument( "First Header entry should be #CHROM: "   + tmp_str ); } break;
            case 1: if ( this->tmp_str != "POS"    ){ throw std::invalid_argument( "Second Header entry should be POS: "     + tmp_str ); } break;
            case 2: if ( this->tmp_str != "ID"     ){ throw std::invalid_argument( "Third Header entry should be ID: "       + tmp_str ); } break;
            case 3: if ( this->tmp_str != "REF"    ){ throw std::invalid_argument( "Fourth Header entry should be REF: "     + tmp_str ); } break; 
            case 4: if ( this->tmp_str != "ALT"    ){ throw std::invalid_argument( "Fifth Header entry should be ALT: "      + tmp_str ); } break;
            case 5: if ( this->tmp_str != "QUAL"   ){ throw std::invalid_argument( "Sixth Header entry should be QUAL: "     + tmp_str ); } break;
            case 6: if ( this->tmp_str != "FILTER" ){ throw std::invalid_argument( "Seventh Header entry should be FILTER: " + tmp_str ); } break;
            case 7: if ( this->tmp_str != "INFO"   ){ throw std::invalid_argument( "Eighth Header entry should be INFO: "    + tmp_str ); } break;
            case 8: if ( this->tmp_str != "FORMAT" ){ throw std::invalid_argument( "Ninth Header entry should be FORMAT: "   + tmp_str ); } break;
        }
        
        if (field_index > 8){  sample_names.push_back( this->tmp_str); }
        
        feild_start = field_end+1;
        field_index++;
    } // End of while loop: field_end < line.size()
    
    //this->nfield_ = field_index;
    
    this->set_nsam( (int)sample_names.size() );
    assert( print_sample_name() );       
}



