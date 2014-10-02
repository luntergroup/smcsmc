#include "segdata.hpp"
using namespace std;

Segment::Segment( string file_name ) { 
    this->file_name_ = file_name;
    //this->FileType = this->file_name_.size() > 0 ? FileType_in : EMPTY; 
    this->init();
}

void Segment::init(){
    in_file.open( this->file_name_.c_str());
    if ( in_file.good() ){
        getline ( in_file, this->tmp_line ); 
        //dout << this->tmp_line <<endl;
        while ( this->tmp_line.size()>0 ){
            buffer_lines.push_back ( this->tmp_line );
            getline ( in_file, this->tmp_line );
        }    
    
    }
    else {
        throw std::invalid_argument ( std::string("Invaild file: ") + this->file_name_ );
        }
    in_file.close();
    this->current_line_index_ = 0;
}

void VariantReader::initialize_read_newLine(){
    allelic_state_at_Segment_start.clear();

    this->feild_start = 0;
    this->field_end   = 0;
    this->field_index = 0;
    this->skip_tmp_line = false;
    
    this->tmp_line = this->buffer_lines[this->current_block_line_];

    }

void Segment::read_new_line(){
    /*! Read VariantReader data, extract mutation site and haplotype
     */ 
    this->current_line_index_++;
    
    //if ( this->FileType == EMPTY ){
        //this->set_empty_line_entry();
        //return;    
    //}
    
    this->initialize_read_newLine(); 

    while ( field_end < this->tmp_line.size() && !this->skip_tmp_line ){
        field_end = min ( this->tmp_line.find('\t',feild_start), this->tmp_line.find('\n', feild_start) );
        this->tmp_str = this->tmp_line.substr( feild_start, field_end - feild_start );
        
        if      ( field_index == 0 ) {    }  
        else if ( field_index == 1 ) { this->extract_field_POS ();    }  /*! Skip the line if two sites are two close, apply the filtering */
        else if ( field_index == 2 ) { this->extract_field_ID ();   }        
        else if ( field_index == 3 ) { this->extract_field_REF () ;   }  /*! REF field, Skip the line, if the length of REF string is greater than one, treated as deletion or indel*/
        else if ( field_index == 4 ) { this->extract_field_ALT () ;   }  /*! ALT field, Skip the line if any alter field */
        else if ( field_index == 5  ) { this->extract_field_VARIANT(); }
        
        feild_start = field_end+1;        
        field_index++;
        }
    }

