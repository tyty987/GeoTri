#ifndef _MACRO_HPP_
#define _MACRO_HPP_

/******************************
* Body
******************************/
#ifdef _TEST_
#include "../include/catch.hpp"
#define _PRIVATE public
#else
#define _PRIVATE public 
#endif

const int MPI_MASTER = 0;


#endif
