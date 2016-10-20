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
        this->src      = "\033[1;31m" + str + "\033[0m";
        this->reason   = "";
    }
    virtual ~InvalidInput() throw() {}
    virtual const char* what () const noexcept {
        return throwMsg.c_str();
    }
};


struct NotEnoughArg : public InvalidInput{
    NotEnoughArg( string str ):InvalidInput( str ){
        this->reason = "Not enough parameters when parsing option: ";
        throwMsg = this->reason + this->src;
    }
    ~NotEnoughArg() throw() {}
};


struct UnknowArg : public InvalidInput{
  UnknowArg( string str ):InvalidInput( str ){
    this->reason = "Unknow option: ";
    throwMsg = this->reason + this->src;
  }
  ~UnknowArg() throw() {}
};

struct FlagsConflict : public InvalidInput{
  FlagsConflict( string str1, string str2 ):InvalidInput( str1 ){
    this->reason = "Flag: ";
    throwMsg = this->reason + this->src + string(" conflict with flag ") + str2;
  }
  ~FlagsConflict() throw() {}
};


struct OutOfEpochRange : public InvalidInput{
  OutOfEpochRange( string str1, string str2 ):InvalidInput( str1 ){
    this->reason = "Problem: epochs specified in -xr/-xc options out of range: ";
    this->src      = "\033[1;31m" + str1 + string(" is greater than ") + str2 + "\033[0m";
    throwMsg = this->reason + src;
  }
  ~OutOfEpochRange() throw() {}
};


struct OutOfRange : public InvalidInput{
  OutOfRange( string str1, string str2 ):InvalidInput( str1 ){
    this->reason = "Flag \"";
    throwMsg = this->reason + this->src + string(" ") + str2 + string("\" out of range [0, 1].");
  }
  ~OutOfRange() throw() {}
};


struct WrongType : public InvalidInput{
  WrongType( string str ):InvalidInput( str ){
    this->reason = "Wrong type for parsing: ";
    throwMsg = this->reason + this->src;
  }
  ~WrongType() throw() {}
};

#endif
