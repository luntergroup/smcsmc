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


bool EvolutionaryEvent::print_event() const {
	EventRecorderdout << "";
#define outstream cout
//#define outstream cout
	if (is_recomb()) {
		outstream << "Event Recomb w=" << weight << " [" << start_base_ << "-" << end_base_ << ")  t=" << setw(5) << start_height << "-" << setw(5) << end_height;
		if (is_recomb_event()) {
			outstream << " (**EVENT**";
			if (event_data == -1) outstream << " xpos=" << a.recomb_pos << ")";
			if (event_data == 0) outstream << " tpos=" << a.recomb_pos << ")";
		}
	} else {
		outstream << "Event CoalMigr w=" << weight << " [" << end_base_ << "] pop=" << a.coal_migr_population << " t=" << setw(5) << start_height << "-" << setw(5) << end_height;
		if (is_coal_event()) 
			outstream << " (**COAL**) ";
		else if (is_migr_event())
			outstream << " (**MIGR** " << event_data << ") ";
	}
	outstream << endl;
    return true;
}


// friend functions

/* Removes event after we're done updating the relevant counters; return true if event was actually deleted */
bool remove_event( EvolutionaryEvent** eventptr_location ) {

	assert (eventptr_location != NULL);
	EvolutionaryEvent* event = *eventptr_location;
	assert (event->children_updated < 0);
	*eventptr_location = event->parent();              // (1) overwrite ptr to event with (2) ptr to parent, so...
	if (event->parent() != NULL) {
		event->parent()->increase_refcount();          // ...(2) increase parent's refcount and
	}
	if (event->decrease_refcount_is_zero()) {          // ...(1) decrease refcount of event
        delete event;                                  // ...which may decrease parent's refcount again (but not to 0)
        return true;                                   // signal: even has been deleted
	}
	return false;									   // signal: even has not been deleted
}

/* Purges previously removed events, and returns first active event (if any) */
EvolutionaryEvent* purge_events( EvolutionaryEvent** eventptr_location ) {

	assert (*eventptr_location != NULL);
	while (1) {
		EvolutionaryEvent* event = *eventptr_location;
		if (event == NULL || !event->is_removed())
			return event;
		remove_event( eventptr_location );
	}
}
