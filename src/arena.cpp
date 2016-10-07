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

#include "arena.hpp"
#include "coalevent.hpp"
#include <iostream>


void Arena::init_arena() {
	
	assert (globalArena == NULL);
	globalArena = this;
	numAllocs = 0;
	numDeallocs = 0;
	maxInUse = 0;
	numSkips = 0 * num_epochs; // to avoid 'variable not used' warning by clang++

	alloc_block();
		
	block_in_use = 0;          // point to block currently being filled
	block_idx = block_size-1;  // first (likely) empty slot in block

}


Arena::~Arena() {

	std::clog << "Arena: allocated " << numAllocs << " blocks, deallocated " << numDeallocs << " blocks, max in use " << maxInUse << " blocks, num skips " << numSkips << std::endl;
	for (int i=0; i<blocks.size(); i++)
		free( allocblocks[i] );

}



void Arena::alloc_block() {

	try {
		void* ptr = malloc( 64 + block_size * sizeof( EvolutionaryEvent ) );
		allocblocks.push_back( ptr );
		ptr = (void*)(64L + ((uintptr_t)ptr & -64L));
		blocks.push_back( ptr );
	} catch ( std::bad_alloc &ba ) {
		std::cerr << "Could not allocate memory (require block of size " << block_size * sizeof( EvolutionaryEvent ) / 1024 << " kb)" << std::endl;
		throw;
	}

	EvolutionaryEvent* block = (EvolutionaryEvent*)(blocks.back());
	for (int j=0; j<block_size; j++) {
		*(double*)(&block[j]) = -1.0;  // mark as 'empty'
	}

}

void Arena::set_next_block() {

	if ( numAllocs - numDeallocs > _arena_max_fill_factor * block_size * blocks.size() ) {
		// blocks are too full -- allocate new
		block_in_use = blocks.size();
		alloc_block();
	} else {
		// move to next block
		block_in_use = (block_in_use + 1) % blocks.size();
	}
}

// static function
void* Arena::allocate( size_t epoch_idx ) {

	Arena* ths = Arena::globalArena;
	ths->numAllocs++;
	ths->maxInUse = std::max( ths->maxInUse, ths->numAllocs - ths->numDeallocs );
	
	EvolutionaryEvent* block = (EvolutionaryEvent*)(ths->blocks[ ths->block_in_use ]);
	int idx = ths->block_idx;
	int block_size = ths->block_size;
	// find first empty spot
	while ( *(double*)(&block[idx]) >= 0.0 ) {
		ths->numSkips++;
		idx--;
		if (idx < 0) {
			ths->set_next_block();
			block = (EvolutionaryEvent*)(ths->blocks[ ths->block_in_use ]);
			idx = block_size - 1;
		}
	}
	ths->block_idx = (idx - 1 + block_size) % block_size;
	if (ths->block_idx == block_size-1)
		ths->set_next_block();
	return (void*)(&block[idx]);	
}

// static function
void Arena::deallocate( void* memptr, size_t epoch_idx ) {

	Arena::globalArena->numDeallocs++;	
	*(double*)memptr = -1.0;   // mark as 'empty'

}
