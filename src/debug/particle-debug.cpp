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


#include "../particleContainer.hpp"

bool ParticleContainer::check_state_orders(){
	dout << "check particle orders, there are " << this->particles.size()<<" particles" <<endl;
	for (size_t i = 0; i < this->particles.size(); i++){
		dout << "Particle " << i << " next genealogy change occurs at position: " << std::setw(14) << this->particles[i]->next_base();
		dout << "  lambda=" << std::setw(10) << particles[i]->getLocalTreeLength();
		dout << " weight =" << std::setw(10) << this->particles[i]->weight() <<std::endl;
        }
    dout << std::endl;
    return true;
    }


void ParticleContainer::print(){
	for (size_t i = 0 ; i < this->particles.size();  i++){
		this->particles[i]->printTree();
		cout << "paritcle " << i<<" has weight " << this->particles[i]->weight()<<endl;
        }
    }


//bool ForestState::print_Coalevent(){
    //cout << " ### Coalescent events:" << endl;
	//for (size_t i = 0 ; i < CoaleventContainer.size() ; i++ ){
		//cout 
        ////<< setw(10) << this->CoaleventContainer[i]->start_height()  << " to " 
             ////<< setw(10) << this->CoaleventContainer[i]->end_height()    << ", " 
             //<< setw(13) << this->CoaleventContainer[i]->opportunity()   << " opportunity for " 
             //<< setw(2)  << this->CoaleventContainer[i]->num_event()     << " coalescent, ";
        //if ( this->CoaleventContainer[i]->event_state() == NOEVENT ){
            //cout<< " potetial coalsecent";
            //}
        //cout << endl;
        //}
	//return true;
    //}
    

//bool ForestState::print_Recombevent(){
    //cout << " ### Recombination events:" << endl;
	//for (size_t i = 0 ; i < RecombeventContainer.size() ; i++ ){
		//cout 
            ////<< setw(10) << this->RecombeventContainer[i]->start_height()  << " to " 
             ////<< setw(10) << this->RecombeventContainer[i]->end_height()    << ", " 
             //<< setw(13) << this->RecombeventContainer[i]->opportunity()   << " opportunity for " 
             //<< setw(2)  << this->RecombeventContainer[i]->num_event()     << " recombination, ";
        //if ( this->RecombeventContainer[i]->event_state() == NOEVENT ){
            //cout<< " potetial recombination";
            //}
        //cout << endl;
        //}
	//return true;
    //}


//bool ForestState::print_Migrevent(){
    //dout << " ### Migration events:" << endl;
	//for (size_t i = 0 ; i < MigreventContainer.size() ; i++ ){
		//dout 
        ////<< setw(10) << this->MigreventContainer[i]->start_height()  << " to " 
             ////<< setw(10) << this->MigreventContainer[i]->end_height()    << ", " 
             //<< setw(13) << this->MigreventContainer[i]->opportunity()   << " opportunity for " 
             //<< setw(2)  << this->MigreventContainer[i]->num_event()     << " migration, ";
        //if ( this->MigreventContainer[i]->event_state() == NOEVENT ){
            //dout << "potetial migration, from pop " << this->MigreventContainer[i]->pop_i() << " to some other population "; 
            //}
        //else {
            //dout << "from pop " << this->MigreventContainer[i]->pop_i() << " to pop " << this->MigreventContainer[i]->mig_pop(); 
            //}
        //dout << endl;
        //}
	//return true;
    //}
