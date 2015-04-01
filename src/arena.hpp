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

#ifndef __ARENA__
#define __ARENA__

#include <iostream>
#include <vector>

using std::vector;

// forward declaration
class Arena;

// maximum occupancy
const double _arena_max_fill_factor = 0.25;

// the global memory arena
class Arena {

public:
	Arena( size_t num_epochs, size_t block_size ) : num_epochs(num_epochs), block_size(block_size) { init_arena(); }
	
	~Arena();
	
	static void* allocate( size_t epoch_idx );
	static void deallocate( void* memptr, size_t epoch_idx );
	
private:

	void init_arena();
	void alloc_block();
	void set_next_block();
	
	static class Arena* globalArena;  // the global memory arena; initialized by init_arena()
	int numAllocs;
	int numDeallocs;
	int maxInUse;
	int numSkips;
	
	vector<void*> allocblocks;
	vector<void*> blocks;
	int block_idx;
	int block_in_use;
	size_t num_epochs;
	size_t block_size;

};

#endif
