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


#include"coalevent.hpp"


bool EvolutionaryEvent::append_event( const EvolutionaryEvent& e )
{
  assert (e.ref_counter == 1);
  if (ref_counter != 1) return false;             // can't add to record that's in multiple use
  if (!is_no_event()) return false;               // this must not have event (can't have 2 events in single object)
  if (is_recomb() != e.is_recomb()) return false; // type of records must match
  if (is_recomb()) {
    if (start_height == e.start_height &&
	    end_height == e.end_height &&
	    end_base_ == e.start_base_ &&
	    weight == e.weight) { // add sequence-wise
      end_base_ = e.end_base_;
      event_data = e.event_data;
      a.recomb_pos = e.a.recomb_pos;
      return true;
    }
    if (start_base_ == e.start_base_ &&
	    end_base_ == e.end_base_ &&
	    end_height == e.start_height &&
	    weight == e.weight) { // add time-wise
      end_height = e.end_height;
      event_data = e.event_data;
      a.recomb_pos = e.a.recomb_pos;
      return true;
    }
  } else { // is_coalmigr
    if (end_base_ == e.end_base_ &&
	    end_height == e.start_height &&
	    a.coal_migr_population == e.a.coal_migr_population &&
	    weight == e.weight) {
      end_height = e.end_height;
      event_data = e.event_data;
      return true;
    }
  }
  return false;
}


bool EvolutionaryEvent::print_event() {
	EventRecorderdout << "";
	if (is_recomb()) {
		dout << "Recombination ";
		if (is_recomb_event()) 
			dout << "(event) ";
	} else {
		dout << "Coal/migration ";
		if (is_coal_event()) 
			dout << "(coal) ";
		else if (is_migr_event())
			dout << "(migr) ";
	}
	dout << setw(10) << this->start_height  << " to ";
    dout << setw(10) << this->end_height    << endl;
    return true;
}


