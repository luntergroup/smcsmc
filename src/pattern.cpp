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
 * Extract number from expression
 */    
size_t extract_Number(const char*& expr) {
    char* end_ptr;
    double res = strtod(expr, &end_ptr); // convert a string to double
    // Advance the pointer and return the result
    expr = end_ptr;
    return (size_t)res;
    }
    

/*!
 * Extract segment factor from expression
 */
size_t extract_SegmentFactors(const char*& expr, vector <size_t> & seg_level1_vec, vector <size_t> & seg_level2_vec) {
    size_t seg_level1 = extract_Number(expr);
    seg_level1_vec.push_back(seg_level1);
    char op = *expr;

    if( op == '+' || op == '\0' ){
        seg_level1_vec.pop_back();
        seg_level1_vec.push_back( (size_t)1 );
        seg_level2_vec.push_back(seg_level1);
        return seg_level1;
        }
    expr++;
    size_t seg_level2 = extract_Number(expr);
    seg_level2_vec.push_back(seg_level2);
    seg_level1 *= seg_level2;
    return seg_level1;
    }


/*!
 * Extract segments from expression
 */
size_t extract_NumberOfSegment (const char*& expr, vector <size_t> & seg_level1_vec, vector <size_t> & seg_level2_vec) {
    size_t num_seg = extract_SegmentFactors(expr, seg_level1_vec, seg_level2_vec);
    for(;;) {
        char op = *expr;
        if( op != '+' ){
            return num_seg;
            }
        expr++;
        size_t next_num_seg = extract_SegmentFactors(expr, seg_level1_vec, seg_level2_vec);        
        num_seg += next_num_seg;
        }
    return num_seg;
    }


/*!
 * PSMC time segment scheme. The time interval [0, top_t] is devided into num_seg number of segments
 */
vector <double> extract_Segment (size_t num_seg, double top_t){
    vector <double> t_i(num_seg);
    for (size_t i = 0; i < num_seg ; i++ ){
        t_i [i] = 0.1 * exp ( (double)i / (double)(num_seg-1) * log(1 + 10 * top_t) ) - 0.1; 
        //cout<<t_i[i]<<" ";
        }
    //cout<<endl;
    return t_i;
    }


vector <double> regroup_Segment (vector <double> old_seg, vector <size_t> & seg_level1_vec, vector <size_t> & seg_level2_vec) {
    vector <double> t_i;
    size_t index = 0;
    for ( size_t seg_i = 0; seg_i < seg_level1_vec.size(); seg_i++){
        for ( size_t i = 0; i < seg_level1_vec[seg_i]; i++){
                t_i.push_back(old_seg[index]);
                index += seg_level2_vec[seg_i];
            }
        }
    dout<<"t_i size = " << t_i.size()<<endl;    
    return t_i;
    }


string convert_pattern (string pattern, double top_t){
    if ( pattern.size() == 0 ){
        return "";
        }
    
    const char * expr = pattern.c_str();
    vector <size_t> seg_level1_vec;
    vector <size_t> seg_level2_vec;
    size_t num_seg = extract_NumberOfSegment ( expr, seg_level1_vec, seg_level2_vec);
    dout << "Number of time segments: " << num_seg <<endl;
    vector <double> t_i = extract_Segment( num_seg, top_t);
    vector <double> new_ti = regroup_Segment (t_i, seg_level1_vec, seg_level2_vec);
    string Ne_array;    
    for (size_t i = 0; i < new_ti.size(); i++ ){
        Ne_array += " -eN " + to_string( new_ti [i] / 2 ) + " 1";
        }
    
    dout<<Ne_array<<endl;
    return Ne_array;    
    }
