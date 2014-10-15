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

//#include "../src/variantReader.hpp"
#include"variantReader.hpp"

void VariantReader::print_vcf_line( ){
    dout << "CHROM = " <<chrom_<<endl;
    dout << "POS = " << site_ <<endl;
    for (size_t i=0;i<sample_names.size();i++){
        dout << sample_names[i]<<": ";
        for (size_t j=0;j<sample_alt[i].size();j++){
            dout << sample_alt[i][j]<<" ";
            }  dout << endl;
        }  dout << endl;
    }


void VariantReader::print(){
    for (size_t i = 0; i < int_vec_of_sample_alt.size(); i++){
            cout << int_vec_of_sample_alt[i]<<" ";
        } cout << endl;
    }
        
    
bool VariantReader::print_sample_name(){
    dout << "Sample names:" << endl;
    for (size_t i = 0; i < this->nsam(); i++){
        dout << sample_names[i]<<" ";
        }  dout << endl;
    return true;
    }
