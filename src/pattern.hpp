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

#include <iostream>
#include <cstring> // strlen
#include <vector>
#include <math.h>
#include <string>


#include <sstream>
#include <iterator>

#ifndef NDEBUG
#define dout std::cout
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#endif

using namespace std;

size_t extract_Number(const char*& expr);
size_t extract_SegmentFactors(const char*& expr, vector <size_t> & seg_level1_vec, vector <size_t> & seg_level2_vec) ;
size_t extract_NumberOfSegment (const char*& expr, vector <size_t> & seg_level1_vec, vector <size_t> & seg_level2_vec) ;
vector <double> extract_Segment (size_t num_seg, double top_t);
vector <double> regroup_Segment (vector <double> old_seg, vector <size_t> & seg_level1_vec, vector <size_t> & seg_level2_vec);
string convert_pattern (string pattern, double top_t);
char *convert(const std::string & s);
