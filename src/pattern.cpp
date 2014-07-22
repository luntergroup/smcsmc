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

#include "pattern.hpp"


/*! 
 * \brief Extract number from expression
 */    
size_t Pattern::extract_Number( ) {
    char* end_ptr;
    double res = strtod(expr_, &end_ptr); // convert a string to double
    // Advance the pointer and return the result
    this->expr_ = end_ptr;
    return (size_t)res;
    }
    

/*!
 * \brief Extract segment factor from expression
 */
size_t Pattern::extract_SegmentFactors( ) {
    size_t seg_level1 = this->extract_Number( );
    this->seg_level1_vec_.push_back(seg_level1);

    this->check_pattern();

    char op = *this->expr_;
    if( op == '+' || op == '\0' ){
        this->seg_level1_vec_.pop_back();
        this->seg_level1_vec_.push_back( (size_t)1 );
        this->seg_level2_vec_.push_back(seg_level1);
        return seg_level1;
        }
    this->expr_++;
    size_t seg_level2 = this->extract_Number( );
    this->seg_level2_vec_.push_back(seg_level2);
    seg_level1 *= seg_level2;
    return seg_level1;
    }


/*!
 * \brief Extract segments from expression
 */
void Pattern::extract_NumberOfSegment ( ) {
    this->num_seg_ = this->extract_SegmentFactors( );
    for(;;) {
        char op = *this->expr_;
                
        if( op != '+' ){
            return;
            }
        this->expr_++;
        this->num_seg_ += this->extract_SegmentFactors( );
        }
    }


/*!
 * \brief PSMC time segment scheme. The time interval [0, top_t] is devided into num_seg number of segments
 */
vector <double> Pattern::extract_Segment ( ){ // This should include the number of samples, as the time intervals should be different according to the number of samples
    vector <double> t_i( this->num_seg_ );
    for (size_t i = 0; i < this->num_seg_ ; i++ ){
        /*! \todo alternative ways to divide the exponential segments */
        // This needs more work, 
        t_i [i] = 0.1 * exp ( (double)i / (double)(this->num_seg_ - 1) * log( 1 + 10 * this->top_t_ ) ) - 0.1; 
        //cout<<t_i[i]<<" ";
        }
    //cout<<endl;
    return t_i;
    }


vector <double> Pattern::regroup_Segment ( vector <double> old_seg ) {
    vector <double> t_i;
    size_t index = 0;
    for ( size_t seg_i = 0; seg_i < this->seg_level1_vec_.size(); seg_i++){
        for ( size_t i = 0; i < this->seg_level1_vec_[seg_i]; i++){
                t_i.push_back(old_seg[index]);
                index += this->seg_level2_vec_[seg_i];
            }
        }
    //dout<<"t_i size = " << t_i.size()<<endl;    
    return t_i;
    }


Pattern::Pattern (string pattern, double top_t):top_t_(top_t){
    if ( pattern[0] == '+' || pattern[0] == '-' ){
        throw std::invalid_argument( string(" Illegal pattern! "));
        }
    
    this->expr_ = pattern.c_str();    
    this->extract_NumberOfSegment ( );
    
    if ( this->num_seg_ < 2 ){
        this->pattern_str = "";
        return;
        }
        
    //dout << "Number of time segments: " << this->num_seg_ <<endl;
    vector <double> t_i = extract_Segment( );
    vector <double> new_ti = regroup_Segment (t_i);

    for (size_t i = 0; i < new_ti.size(); i++ ){
        this->pattern_str += " -eN " + to_string( new_ti [i] / 2 ) + " 1";
        }
    
    //dout << this->pattern_str << endl;    
    }


void Pattern::check_pattern ( ){
    char op = *this->expr_;
    if ( op != '+' && op != '\0' && op != '*' && !isdigit(op) ){
        throw std::invalid_argument( string("Character ") + op + string(" is invalid for pattern!"));
        }    
    char op_next = *(this->expr_ + 1);    
    if ( op == '+' && ( op_next == '*' || op_next == '+' || op_next == '\0' ) ){
        throw std::invalid_argument( string(" Illegal pattern! "));
        }
    if ( op == '*' && ( op_next == '*' || op_next == '+' || op_next == '\0' || isdigit(op_next) ) ){
        throw std::invalid_argument( string(" Illegal pattern! "));
        }
    }
