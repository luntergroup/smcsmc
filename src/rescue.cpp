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


#include "rescue.hpp"
#include <fstream>      // std::ifstream

RescueHist::RescueHist ( string HIST_filename ) {
    // open file, loop through, if there is "====", rescure again, and add EMstep by one
    ifstream in_file( HIST_filename.c_str() );

    //if (hist){
        //Ne_file << "=========\n"; 
        //}
    //Ne_file << "RE\t" << this->model->recombination_rate() << "\n";
    
    in_file.close();
    }


void RescueHist::rescue_RE () { // this need to rescale by the sequence length ...
    
    }


void RescueHist::rescue_ME () {
    
    }

void RescueHist::rescue_NE () {
    
    }
