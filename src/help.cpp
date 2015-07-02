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

#include "help.hpp"

void Help_option(){
    //cout << "Too few command line arguments" << endl;
    cout << "Options:" << endl;
    cout << setw(10)<<"-Np"     << setw(5) << "INT" << "  --  " << "Number of particles [ 1000 ]" << endl;
    cout << setw(10)<<"-ESS"    << setw(5) << "FLT" << "  --  " << "Fractional ESS threshold for resampling (1 = use random likelihoods) [ 0.6 ]" << endl;
    cout << setw(10)<<"-p"      << setw(5) << "STR" << "  --  " << "Pattern of time segments [ \"3*1+2*3+4\" ]" <<endl;
    cout << setw(10)<<"-tmax"   << setw(5) << "FLT" << "  --  " << "Maximum time, in unit of 4N0 [ 3 ]" <<endl;
    cout << setw(10)<<"-EM"     << setw(5) << "INT" << "  --  " << "EM iterations [ 20 ]" << endl;
    //cout << setw(10)<<"-lag"    << setw(5) << "FLT" << "  --  " << "Lagging step [ 1000 ]." << endl;
    cout << setw(10)<<"-seg"    << setw(5) << "STR" << "  --  " << "Data file in seg format [ Chrom1.seg ]" << endl;
    //cout << setw(20)<<"-buff BUFFSIZE" << "  --  " << "User define the size of buffer for the seg file BUFFSIZE." << endl;
    cout << setw(10)<<"-o"      << setw(5) << "STR" << "  --  " << "Prefix for output files" << endl;
    cout << setw(10)<<"-online" << setw(5) << " "   << "  --  " << "Perform online EM" << endl;
    cout << setw(10)<<"-xr"     << setw(5) << "INT" << "  --  " << "Epoch or epoch range to exclude from recombination EM (1-based, closed)" << endl;
    cout << setw(10)<<"-xc"     << setw(5) << "INT" << "  --  " << "Epoch or epoch range (e.g. 1-10) to exclude from coalescent/migration EM" << endl;
    cout << setw(10)<<"-log"    << setw(5) << " "   << "  --  " << "Generate *.log file" << endl;
    cout << setw(10)<<"-heat"   << setw(5) << " "   << "  --  " << "Generate *TMRCA and *WEIGHT for heatmap" << endl;
    cout << setw(10)<<"-v"      << setw(5) << " "   << "  --  " << "Display timestamp and git versions" << endl;
};


void Help_example(){
    cout << "    Examples:" << endl;
    cout << "smcsmc 10 -nsam 3" << endl;
    cout << "./smcsmc -Np 5 -t 0.002 -r 400 -npop 20000 -seg eg_seg.seg -buff 4" << endl;
    cout << "./smcsmc -Np 5 -t 0.002 -r 400 -npop 20000 -seg eg_seg.seg" << endl;
    cout << "./smcsmc -Np 6 -t 0.0002 -r 30 -npop 10000 -seed 1314 -seg eg_seg.seg" << endl;
    cout << "./smcsmc -Np 7 -t 0.002 -log -r 400 -seg eg_seg.seg " << endl;
};

void Help_header(){
    cout << "smcsmc" << endl;
    cout << "  version: " << VERSION << endl;
    cout << "  authored by Sha (Joe) Zhu and Gerton Lunter " <<endl;
    Help_option();
    Help_example();
    exit(0);
};

void Help_version(string date, string smcsmc, string scrm, std::ostream &output){
    output << "Program was compiled on: " << date << endl;
    output << "smcsmc version: " << smcsmc << endl;
    output << "scrm version:   " << scrm   << endl;
}
