/*
 * smcsmc is short for particle filters for ancestral recombination graphs.
 * This is a free software for demographic inference from genome data with particle filters.
 *
 * Copyright (C) 2013-2017 Donna Henderson, Sha (Joe) Zhu and Gerton Lunter
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

#include <string>  /* string */
#include <exception>

#ifndef EXCEPTION
#define EXCEPTION

using namespace std;

struct InvalidInput : std::exception {
    string src;
    string reason;
    string throwMsg;

    InvalidInput( ){
        this->src      = "";
        this->reason   = "";
    }

    explicit InvalidInput( string str ){
        this->src      = str;
        this->reason   = "";
    }
    virtual ~InvalidInput() throw() {}
    virtual const char* what () const noexcept {
        if (throwMsg.size() == 0) 
            const_cast<InvalidInput*>(this)->throwMsg = src + " " + reason;
        return throwMsg.c_str();
    }
};



#endif
