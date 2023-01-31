#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <iostream>
#include <iomanip>

#include "defs.hpp"

template<typename T> void PrintElement(std::ostream& os, T t, const int& width)
{
  os << std::left << std::setw(width) << std::setfill(' ') << t;
}

void ChangeRunDir(const char *pdir);
Real* LogSpace(Real a, Real b, int k);
void WriteData(char *filename,
               Real *nH, Real *press,
               Real *x_h2, Real *x_hii, Real *x_e,
               Real *chi_pe, Real *chi_lw,
               Real *chi_h2, Real *chi_ci,
               Real *xi_cr,
               const int n);


#endif // UTILS_UTILS_HPP_
