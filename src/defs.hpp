//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
#ifndef DEFS_HPP_
#define DEFS_HPP_

// try/throw/catch C++ exception handling (ENABLE_EXCEPTIONS or DISABLE_EXCEPTIONS)
// (enabled by default)
#define ENABLE_EXCEPTIONS

// use single precision floating-point values (binary32)? default=0 (false; use binary64)
#define SINGLE_PRECISION_ENABLED 0

// primitive type alias that allows code to run with either floats or doubles
#if SINGLE_PRECISION_ENABLED
using Real = float;
#else
using Real = double;
#endif

#define SQR(x) ( (x)*(x) )
#define CUBE(x) ( (x)*(x)*(x) )
#define FOURTH(x) ( (x)*(x)*(x)*(x) )

#define PI 3.1415926535897932
#define ONE_3RD  0.3333333333333333
#define TWO_3RDS 0.6666666666666667
#define FOUR_3RDS 1.333333333333333

#define TINY_NUMBER 1.0e-20
#define HUGE_NUMBER 1.0e+20

#define XEL_MAX 1.2006199779862501

#ifdef ENABLE_EXCEPTIONS
#define ATHENA_ERROR(x) throw std::runtime_error(x.str().c_str())
#else
#define ATHENA_ERROR(x) std::cout << x.str(); std::exit(EXIT_FAILURE)
#endif

#endif // DEFS_HPP_
