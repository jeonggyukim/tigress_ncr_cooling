
// POSIX C extensions
#include <sys/stat.h>  // mkdir()
#include <unistd.h>    // chdir()

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

#include "defs.hpp"
#include "units.hpp"
#include "utils.hpp"

//----------------------------------------------------------------------------------------
//! \fn void ChangeRunDir(const char *pdir)
//! \brief change to input run directory; create if it does not exist yet
void ChangeRunDir(const char *pdir) {
  std::stringstream msg;

  if (pdir == nullptr || *pdir == '\0') return;

  mkdir(pdir, 0775);
  if (chdir(pdir)) {
    msg << "### FATAL ERROR in function [ChangeToRunDir]" << std::endl
        << "Cannot cd to directory '" << pdir << "'";
    ATHENA_ERROR(msg);
  }

  return;
}

Real* LogSpace(Real a, Real b, int k) {
  // std::vector<Real> logspace;
  // logspace.reserve(k);
  Real *logspace = new Real[k]();
  for (int i=0; i<k; ++i) {
    logspace[i] = std::pow(10, a + i*(b-a)/(k-1));
  }
  return logspace;
}



void WriteData(std::string filename,
               Real *nH, Real *press,
               Real *x_h2, Real *x_hii, Real *x_e,
               Real *chi_pe, Real *chi_lw,
               Real *chi_h2, Real *chi_ci,
               Real *xi_cr,
               const int n) {

  const int prec = 5;
  const int w = prec+7;
  Units *punit = new Units();
  
  std::ofstream os;
  os.open(filename);
  os << std::setw(w) << "nH "
     << std::setw(w) << "pok "
     << std::setw(w) << "xH2 "
     << std::setw(w) << "xHII "
     << std::setw(w) << "xe "
     << std::setw(w) << "chi_pe "
     << std::setw(w) << "chi_lw "
     << std::setw(w) << "chi_h2 "
     << std::setw(w) << "chi_ci "
     << std::setw(w) << "xi_cr "
     << std::endl;
  
  // std::cout << std::setprecision(5) << std::left;
  
  os << std::scientific << std::setprecision(prec) << std::left;
  for (int i=0; i<n; ++i) {
    PrintElement(os, nH[i], w);
    PrintElement(os, press[i]*punit->Pressure/Constants::kB, w);
    PrintElement(os, x_h2[i], w);
    PrintElement(os, x_hii[i], w);
    PrintElement(os, x_e[i], w);
    PrintElement(os, chi_pe[i], w);
    PrintElement(os, chi_lw[i], w);
    PrintElement(os, chi_h2[i], w);
    PrintElement(os, chi_ci[i], w);
    PrintElement(os, xi_cr[i], w);
    os << std::endl;
  }

  delete punit;
}
