/** \file
 *
 *  This file contains declarations of a small collection of 
 *  general-purpose functions that are used throughout the npf code.
 */

#ifndef H_BASICS_H
#define H_BASICS_H

#include <iostream>
#include <vector>
#include <cstdlib>

/** This function can be called when a fatal error occurs.
 *  It writes the error message \a msg to the standard error stream,
 *  then terminates the program, returning the value 1 to the 
 *  operating system.
 */
inline void npfFatalError( const char *msg )
{
	std::cerr << msg << std::endl;
	std::exit(1);
}

/** Believe it or not, but C++ does not have a function sqr
 *  for calculating the square of an argument x. Here we provide
 *  it as a so-called function template. As a result this function
 *  can be used for all types that support multiplication: int,
 *  double, complex, ...
 *
 */
template <typename T>
inline T sqr(const T& value)
{
	return value*value;
}

#endif // H_BASICS_H

