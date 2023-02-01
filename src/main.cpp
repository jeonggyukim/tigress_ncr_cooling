//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
#include <cmath>
#include <string>

#include "defs.hpp"
#include "cooling_tigress.hpp"
#include "units.hpp"
#include "utils.hpp"
#include "parameter_input.hpp"

static const Real muH = 1.4;

int main(int argc, char *argv[]) {

#ifdef ATTACH_DEBUGGER
  int ii = 0;
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  std::cout << "PID " << getpid() << " on " << hostname
            << " ready for debugger attach" << std::endl;
  fflush(stdout);
  while (ii == 0)
    sleep(5);
#endif

  char *input_filename = nullptr;
  char *prundir = nullptr;
  int iarg_flag = 0;  // set to 1 if -i <file> argument is on cmdline
  int narg_flag = 0;  // set to 1 if -n        argument is on cmdline
  
  // Check for command line options and respond.
  for (int i=1; i<argc; i++) {
    // If argv[i] is a 2 character string of the form "-?" then:
    if (*argv[i] == '-'  && *(argv[i]+1) != '\0' && *(argv[i]+2) == '\0') {
      // check validity of command line options + arguments:
      char opt_letter = *(argv[i]+1);
      switch(opt_letter) {
        // options that do not take arguments:
        case 'n':
        case 'c':
        case 'h':
          break;
          // options that require arguments:
        default:
          if ((i+1 >= argc) // flag is at the end of the command line options
              || (*argv[i+1] == '-') ) { // flag is followed by another flag
              std::cout << "### FATAL ERROR in main" << std::endl
                        << "-" << opt_letter << " must be followed by a valid argument\n";
              return 0;
          }
      }
      switch(*(argv[i]+1)) {
      case 'i':                      // -i <input_filename>
        input_filename = argv[++i];
        iarg_flag = 1;
        break;
      case 'd':                      // -d <run_directory>
        prundir = argv[++i];
        break;
      case 'n':
        narg_flag = 1;
        break;
      case 'h':
      default:
        std::cout << "tigress_ncr_cooling " << std::endl;
        std::cout << "Usage: " << argv[0] << " [options] [block/par=value ...]\n";
        std::cout << "Options:" << std::endl;
        std::cout << "  -i <file>       specify input file [athinput]\n";
        std::cout << "  -d <directory>  specify run dir [current dir]\n";
        
        return(0);
        break;
      }
    } // else if argv[i] not of form "-?" ignore it here (tested in ModifyFromCmdline)
  }

  if (input_filename == nullptr) {
    // no input file is given
    std::cout << "### FATAL ERROR in main" << std::endl
              << "No input file or restart file is specified." << std::endl;
    return(0);
  }
  ParameterInput *pinput;
  IOWrapper infile;

#ifdef ENABLE_EXCEPTIONS
  try {
#endif
    pinput = new ParameterInput;
    if (iarg_flag == 1) {
      infile.Open(input_filename, IOWrapper::FileMode::read);
      pinput->LoadFromFile(infile);
      infile.Close();
    }
    pinput->ModifyFromCmdline(argc ,argv);
#ifdef ENABLE_EXCEPTIONS
  }
  catch(std::bad_alloc& ba) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "memory allocation failed initializing class ParameterInput: "
              << ba.what() << std::endl;
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message
    return(0);
  }
#endif // ENABLE_EXCEPTIONS

  // Dump input parameters and quit if code was run with -n option.
  if (narg_flag) {
    pinput->ParameterDump(std::cout);
    return(0);
  }
  
  Units *punit = new Units();
  punit->PrintCodeUnits();

  int ngrid = pinput->GetReal("problem", "ngrid");
  Real log_nH_min = pinput->GetReal("problem","log_nH_min");
  Real log_nH_max = pinput->GetReal("problem","log_nH_max");
  int flag_dust_cool = pinput->GetInteger("cooling","flag_dust_cool");
  
  Real z_g = pinput->GetReal("problem","z_g");
  Real z_d = pinput->GetReal("problem","z_d");
  Real chi0 = pinput->GetReal("problem","chi0");
  Real xi_cr0 = pinput->GetReal("problem","xi_cr0");
  
  Real sigma_dust_pe0 = pinput->GetReal("opacity","sigma_dust_pe0");
  Real sigma_dust_lw0 = pinput->GetReal("opacity","sigma_dust_lw0");
  Real sigma_dust_pe_cgs = sigma_dust_pe0*z_d;
  Real sigma_dust_lw_cgs = sigma_dust_lw0*z_d;

  Real nH0 = pinput->GetReal("shielding","nH0");
  Real len_shld0 = pinput->GetReal("shielding","len_shld0");;
  Real pow_idx_shld = pinput->GetReal("shielding","pow_idx_shld");
    
  Real temp_mu_init = pinput->GetReal("problem","temp_mu_init");
  Real x_h2__;
  Real x_h2_, x_hii_, x_hi_, x_e_, x_cii_, x_ci_;
  
   // Shielding column and threshold for cr attenuation
  Real col_shld_cgs;
  Real col_shld_cr_cgs = pinput->GetReal("shielding","col_shld_cr_cgs");

  Real t_end_Myr = pinput->GetReal("problem","t_end_Myr");
  Real t_end = t_end_Myr * Constants::Myr / punit->Time;
  std::cout << "t_end " << t_end << std::endl;
  
  Real *nH = LogSpace(log_nH_min, log_nH_max, ngrid);
  Real *press = new Real[ngrid]();
  Real *x_h2 = new Real[ngrid]();
  Real *x_hii = new Real[ngrid]();
  Real *x_e = new Real[ngrid]();
  Real *chi_pe =  new Real[ngrid]();
  Real *chi_lw =  new Real[ngrid]();
  Real *chi_h2 =  new Real[ngrid]();
  Real *chi_ci =  new Real[ngrid]();
  Real *xi_cr =  new Real[ngrid]();
  Real *len_shld =  new Real[ngrid]();
  
  // int w = 10;
  // std::cout << "|" << std::setw(w) << "nH"
  //           << "|" << std::setw(w) << "T1"
  //           << "|" << std::setw(w) << "pok"
  //           << "|" << std::setw(w) << "xH2"
  //           << "|" << std::setw(w) << "xHII"
  //           << "|" << std::setw(w) << "xe"
  //           << "|" << std::endl;
  // std::cout << std::setprecision(5) << std::left;
  
  // Set initial data
  for (int i=0; i<ngrid; ++i) {
    len_shld[i] = len_shld0*std::pow(nH[i]/nH0, pow_idx_shld);
    col_shld_cgs = nH[i]*len_shld[i]*punit->Length;
    chi_pe[i] = chi0*exp(-col_shld_cgs*sigma_dust_pe_cgs);
    chi_lw[i] = chi0*exp(-col_shld_cgs*sigma_dust_lw_cgs);
    chi_h2[i] = chi0*exp(-col_shld_cgs*sigma_dust_lw_cgs);
    chi_ci[i] = chi0*exp(-col_shld_cgs*sigma_dust_lw_cgs);

    xi_cr[i] = xi_cr0;
    if (col_shld_cgs > col_shld_cr_cgs) {
      xi_cr[i] *= col_shld_cr_cgs/col_shld_cgs;
    }

    // Equilibrium H2 abundance without self-shielding
    x_h2__ = get_xH2(nH[i], muH*temp_mu_init/1.1, z_d, xi_cr[i], chi_h2[i]);
    x_e_ = get_xe_fast(nH[i], muH*temp_mu_init/(1.1 - x_h2__),
                       z_d, z_g, xi_cr[i], chi_pe[i], chi_ci[i], chi_h2[i],
                       x_h2__, 1, &x_hii_, &x_h2_, &x_cii_, &x_ci_);
    x_hi_ = std::max(0.0, 1.0 - 2.0*x_h2_ - x_hii_);
    
    press[i] = (1.1 + x_e_ - x_h2_)*nH[i]*Constants::kB*temp_mu_init/punit->Pressure;
    x_h2[i] = x_h2_;
    x_hii[i] = x_hii_;
    x_e[i] = x_e_;

    // std::cout << " "
    //           << nH[i] << " "
    //           << temp_mu << " "
    //           << press[i]*punit->Pressure/Constants::kB << " "
    //           << x_h2[i] << " "
    //           << x_hii[i] << " "
    //           << x_e[i] << " "
    //           << std::endl;
  }

  CoolingSolverTigress *pcool = new CoolingSolverTigress(flag_dust_cool,
                                                         sigma_dust_pe0,
                                                         sigma_dust_lw0,
                                                         punit);

  
  pcool->OperatorSplitSolver1D(nH, press,
                               x_h2, x_hii, x_e,
                               chi_pe, chi_lw,
                               chi_h2, chi_ci,
                               xi_cr, len_shld,
                               chi0,
                               z_g, z_d,
                               0, ngrid, t_end);

  // for (int i=0; i<n; ++i) {
  //   std::cout << " "
  //             << nH[i] << " "
  //             << press[i]/nH[i]*punit->Temperature_mu << " "
  //             << press[i]*punit->Pressure/Constants::kB << " "
  //             << x_h2[i] << " "
  //             << x_hii[i] << " "
  //             << x_e[i] << " "
  //             << std::endl;
  // }

  ChangeRunDir(prundir);
  std::string filename = pinput->GetString("job","problem_id");
  filename.append(".txt");

  WriteData(filename,
            nH, press,
            x_h2, x_hii, x_e,
            chi_pe, chi_lw, chi_h2, chi_ci,
            xi_cr,
            ngrid);

  delete[] nH;
  delete[] press;
  delete[] x_h2;
  delete[] x_hii;
  delete[] x_e;
  
  delete punit;
  delete pcool;
  
  return(0);
  
}
