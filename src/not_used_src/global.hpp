//global.hpp
//#include<iostream>
#include<string>
#include<sstream>
#include<fstream>
//#include<vector.h>
#include<iomanip>
#include"particle.hpp"
#include<valarray>
#include <boost/assign/std/vector.hpp>

#ifndef NDEBUG
#define dout std::cout
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#endif

