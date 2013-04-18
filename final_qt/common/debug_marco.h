/** debug_macro.h
 ** Brief: This is the header file of all the useful macros for debugging.
 ** Project: large-scale fluids
 ** Date: 04/16/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#ifndef DEBUG_MARCO_H
#define DEBUG_MARCO_H

#include <iostream>
#include <assert.h>


#define dout_var(var) dout << #var" = " << var << endl;
#define dout std::cout << __FILE__ << "@" << __LINE__ << " (" << __PRETTY_FUNCTION__ << "):

#endif // DEBUG_MARCO_H
