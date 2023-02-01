//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cooling_function_tigress.cpp
//! \brief cooling solver for TIGRESS-NCR cooling

// C headers

// C++ headers
#include <iostream>
#include <cmath>
#include <algorithm>

#include "units.hpp"
#include "parameter_input.hpp"
#include "cooling_tigress.hpp"

static const int flag_cool_hyd_cie = 0;

// static const Real gamma = 5.0/3.0;
// static const Real gm1 = 2.0/3.0;
// static const Real igm1 = 3.0/2.0;
static const Real temp_mu_floor = 2.0;

static const Real x_h2_cut = 0.0;
static const Real x_hi_cut = 0.0;
static const Real x_hii_cut = 0.0;

static const int flag_h2_cooling = 1;
static const int flag_h2_heating = 1;

static const Real b5_inv=0.3333333333333333;

static const Real u_rad_pe_isrf_cgs = 7.613e-14;
static const Real u_rad_lw_isrf_cgs = 1.335e-14;
static const Real xi_diss_h2_isrf = 5.7e-11;

static const Real xCstd = 1.6e-4;
static const Real xOstd = 3.2e-4;

static const Real temp_dust0 = 5.0;

static const Real dvdr = 3.240779289444365e-14;

CoolingSolverTigress::CoolingSolverTigress(ParameterInput *pin,
                                           Units *punit) :
  punit(punit),
  flag_cool_dust_(pin->GetInteger("cooling","flag_cool_dust")),
  flag_cool_hyd_cie_(pin->GetOrAddInteger("cooling","flag_cool_hyd_cie",0)),
  nsub_max_(pin->GetInteger("cooling","nsub_max")),
  cfl_cool_sub_(pin->GetOrAddReal("cooling","cfl_cool_sub",0.1)),
  temp_hot0_(pin->GetReal("cooling","temp_hot0")),
  temp_hot1_(pin->GetReal("cooling","temp_hot1")),
  sigma_dust_pe0_(pin->GetReal("opacity","sigma_dust_pe0")),
  sigma_dust_lw0_(pin->GetReal("opacity","sigma_dust_lw0")),
  muH_(1.4),
  code_den_to_nH_(punit->Density/(muH_*Constants::mH)),
  gm1_(pin->GetOrAddReal("problem","gamma",5.0/3.0) - 1.0),
  igm1_(1.0/gm1_)
{
}

CoolingSolverTigress::~CoolingSolverTigress() {};

void CoolingSolverTigress::OperatorSplitSolver1D(Real *den, Real *press,
                                                 Real *x_h2, Real *x_hii, Real *x_e,
                                                 Real *chi_pe, Real *chi_lw,
                                                 Real *chi_h2, Real *chi_ci,
                                                 Real *xi_cr, Real *len_shld,
                                                 const Real chi0,
                                                 const Real z_g, const Real z_d,
                                                 int il, int iu, const Real dt) {
  
  // Real dt = pmb->pmy_mesh->dt;

#pragma omp simd
  for (int i=il; i<iu; ++i) {
    CoolVar cv = {};
    SetupVariables1D(den, press,
                     x_h2, x_hii, x_e,
                     chi_pe, chi_lw,
                     chi_h2, chi_ci,
                     xi_cr, len_shld,
                     chi0, z_g, z_d,
                     cv, i);
    CoolingExplicitSubcycling(cv, dt);
    UpdateVariables1D(den, press,
                      x_h2, x_hii, x_e,
                      chi_h2, chi_ci,
                      cv, i);
  }
  
  return;
}


// void CoolingSolverTigress::OperatorSplitSolver(MeshBlock *pmb) {
// // void CoolingSolverTigress::OperatorSplitSolver() {

//   // Set index bounds. Since boundary communication will not be called, need to
//   // solve cooling in the ghost zones
//   int il = pmb->is - NGHOST;
//   int iu = pmb->ie + NGHOST;
//   int jl = pmb->js;
//   int ju = pmb->je;
//   if (pmb->pmy_mesh->f2) {
//     jl -= (NGHOST);
//     ju += (NGHOST);
//   }
//   int kl = pmb->ks;
//   int ku = pmb->ke;
//   if (pmb->pmy_mesh->f3) {
//     kl -= (NGHOST);
//     ku += (NGHOST);
//   }

//   // int il = 0;
//   // int iu = 10;
//   // int jl = 0;
//   // int ju = 0;
//   // int kl = 0;
//   // int ku = 0;
//   // Real dt = 0.1;
  
//   CoolVar cv;


//   Real dt = pmb->pmy_mesh->dt;

//   for (int k=kl; k<=ku; ++k) {
//     for (int j=jl; j<=ju; ++j) {
// #pragma omp simd
//       for (int i=il; i<=iu; ++i) {
//         SetupVariables(pmb, cv, k, j, i);
//         CoolingExplicitSubcycling(cv, dt);
//         UpdateHydroVariables(pmb, cv, k, j, i);        
//       }
//     }
//   }

// }

void CoolingSolverTigress::SetupVariables1D(Real *den, Real *press,
                                            Real *x_h2, Real *x_hii, Real *x_e,
                                            Real *chi_pe, Real *chi_lw,
                                            Real *chi_h2, Real *chi_ci,
                                            Real *xi_cr, Real *len_shld,
                                            const Real chi0,
                                            const Real z_g, const Real z_d,
                                            CoolVar& cv,
                                            const int i) {
  
  Real xi_ph_hi = 0.0;
  Real xi_ph_h2 = 0.0;
  
  // metallicity and dust abundance
  cv.z_g = z_g;
  cv.z_d = z_d;

  // Real w_d  = den[i];
  // Real w_p  = press[i];
  // Real u_e  = u_e[i];

  Real x_h2_ = x_h2[i];
  Real x_hii_ = x_hii[i];
  Real x_e_ = x_e[i];
  Real x_e_max_ = 1.20061997798625 + 0.00961998*(cv.z_g - 1.0);
  Real x_hi_;
    
  // Store density, non-thermal energy, and initial pressure
  cv.den = den[i];
  cv.nH = cv.den*code_den_to_nH_;
  cv.press_before = press[i];

  // Store pressure floor and non-thermal (e.g., kinetic and magnetic) energy
  cv.press_floor = cv.den*temp_mu_floor/punit->Temperature_mu;
  // cv.e_non_thermal = u_e - w_p*igm1_;

  // Set hydrogen abundances
  // Enforce closure relation (similar to approach taken by RAMSES)
  Real x_sum_ = x_hii_ + 2.0*x_h2_;
  if (x_sum_ > 1.0) {
    x_hii_ = x_hii_/x_sum_;
    x_h2_ = x_h2_/x_sum_;
  }
  if (x_h2_ < 0.25) {
    x_h2_ = std::max(TINY_NUMBER, x_h2_);
    x_h2_ = std::min(0.5, x_h2_);
    x_hii_ = std::max(TINY_NUMBER, x_hii_);
    x_hii_ = std::min(1.0 - 2.0*x_h2_, x_hii_);
  } else {
    x_hii_ = std::max(TINY_NUMBER, x_hii_);
    x_hii_ = std::min(1.0, x_hii_);
    x_h2_ = std::max(TINY_NUMBER, x_h2_);
    x_h2_ = std::min(1.0 - x_hii_, x_h2_);
  }
  x_hi_ = std::max(TINY_NUMBER, 1.0 - x_hii_ - 2.0*x_h2_);
  
  x_e_ = std::max(TINY_NUMBER, x_e_);
  x_e_ = std::min(x_e_max_, x_e_);

  cv.x_h2 = x_h2_;
  cv.x_hii = x_hii_;
  cv.x_hi = x_hi_;
  cv.x_e = x_e_;

  // Apply pressure floor before solving the cooling
  cv.dpress_heat = 0.0;
  cv.dpress_cool = 0.0;
  cv.dpress_floor = 0.0;
  if (cv.press_before < cv.press_floor) {
    // store heating added by flooring
    cv.dpress_floor += cv.press_floor - cv.press_before;
  }
  cv.press_before = std::max(cv.press_before, cv.press_floor);
  cv.press = cv.press_before;

  // Set temp_mu and temp
  cv.temp_mu = cv.press/cv.den*punit->Temperature_mu;
  cv.mu = muH_/(1.1 + cv.x_e - cv.x_h2);
  cv.temp  = cv.mu*cv.temp_mu;
  cv.temp_ = cv.temp*(1.0 + dlntemp_);
  
  // SetRadiationVariables(pmb, cv);

  // Note that chi_fuv is the normalized intensity in the FUV (PE+LW) band used
  // for PE heating and chi_pe is the normalized intensity in the PE band.
  cv.chi_fuv = (chi_pe[i]*u_rad_pe_isrf_cgs + chi_lw[i]*u_rad_lw_isrf_cgs)/
    (u_rad_pe_isrf_cgs + u_rad_lw_isrf_cgs);
  cv.chi_lw = chi_lw[i];
  cv.chi_h2 = chi_h2[i];
  cv.chi_ci = chi_ci[i];
  cv.chi_co = chi_lw[i];
  
  // SetCosmicRayVariables(pmb, cv);
  cv.xi_cr = xi_cr[i];
  cv.len_shld = len_shld[i];

  cv.xi_ph_hi = xi_ph_hi;
  cv.xi_ph_h2 = xi_ph_h2;
  cv.xi_diss_h2 = chi_h2[i]*xi_diss_h2_isrf;
  
  cv.dvdr = dvdr;

  if (flag_cool_dust_) {
    // TODO add LyC
    cv.heating_dust_uv =
      (cv.z_g*sigma_dust_pe0_)*(chi_pe[i]*u_rad_pe_isrf_cgs) + 
      (cv.z_g*sigma_dust_lw0_)*(chi_lw[i]*u_rad_lw_isrf_cgs);
    // Additional heating by other (e.g., optical) photons
    // CAUTION: minimum heating is assumed to be ISRF (not proportional to chi0)
    // this assumption was adopted in Athena-TIGRESS code
    Real tau_pe_eff = -std::log(chi_pe[i]/chi0);
    Real heating_dust_min =
      ((cv.z_g*sigma_dust_pe0_)*u_rad_pe_isrf_cgs +
       (cv.z_g*sigma_dust_lw0_)*u_rad_lw_isrf_cgs)*exp(-tau_pe_eff/1.87);
    cv.heating_dust_uv += heating_dust_min;
    cv.heating_dust_uv *= Constants::c;
    cv.temp_dust = temp_dust0;
  }
  
  return;
}


// void CoolingSolverTigress::SetupVariables(MeshBlock *pmb, CoolVar& cv,
//                                                  const int k, const int j, const int i) {

//   Real z_g = 1.0;
//   Real z_d = 1.0;
  
//   Real chi_pe = 1.0;
//   Real chi_ci = 1.0;
//   Real chi_co = 1.0;
//   Real xi_diss_h2 = 0.0;
//   Real xi_ph_hi = 0.0;
//   Real xi_ph_h2 = 0.0;
//   Real xi_cr = 2.0e-16;
  
//   Real w_d  = pmb->phydro->w(IDN,k,j,i);
//   Real w_p  = pmb->phydro->w(IPR,k,j,i);
//   Real u_e  = pmb->phydro->u(IEN,k,j,i);

//   Real x_h2 = pmb->pscalars->r(0,k,j,i);
//   Real x_hii = pmb->pscalars->r(1,k,j,i);
//   Real x_e = pmb->pscalars->r(2,k,j,i);
//   Real x_e_max = 1.20061997798625 + 0.00961998*(z_g - 1.0);
//   Real x_hi;
  
//   // metallicity and dust abundance
//   cv.z_g = z_g;
//   cv.z_d = z_d;
  
//   // Store density, non-thermal energy, and initial pressure
//   cv.den = w_d;
//   cv.nH = cv.den*code_den_to_nH_;
//   cv.press_before = w_p;

//   // Store pressure floor and non-thermal (e.g., kinetic and magnetic) energy
//   cv.press_floor = cv.den*temp_mu_floor/punit->Temperature_mu;
//   cv.e_non_thermal = u_e - w_p*igm1_;

//   // Set hydrogen abundances
//   // Enforce closure relation (similar to approach taken by RAMSES)
//   Real x_sum = x_hii + 2.0*x_h2;
//   if (x_sum > 1.0) {
//     x_hii = x_hii/x_sum;
//     x_h2 = x_h2/x_sum;
//   }
//   if (x_h2 < 0.25) {
//     x_h2 = std::max(TINY_NUMBER, x_h2);
//     x_h2 = std::min(0.5, x_h2);
//     x_hii = std::max(TINY_NUMBER, x_hii);
//     x_hii = std::min(1.0 - 2.0*x_h2, x_hii);
//   } else {
//     x_hii = std::max(TINY_NUMBER, x_hii);
//     x_hii = std::min(1.0, x_hii);
//     x_h2 = std::max(TINY_NUMBER, x_h2);
//     x_h2 = std::min(1.0 - x_hii, x_h2);
//   }
//   x_hi = std::max(TINY_NUMBER, 1.0 - x_hii - 2.0*x_h2);
  
//   x_e = std::max(TINY_NUMBER, x_e);
//   x_e = std::min(x_e_max, x_e);

//   cv.x_h2 = x_h2;
//   cv.x_hii = x_hii;
//   cv.x_hi = x_hi;
//   cv.x_e = x_e;

//   // Apply pressure floor before solving the cooling
//   cv.dpress_heat = 0.0;
//   cv.dpress_cool = 0.0;
//   cv.dpress_floor = 0.0;
//   if (cv.press_before < cv.press_floor) {
//     // store heating added by flooring
//     cv.dpress_floor += cv.press_floor - cv.press_before;
//   }
//   cv.press_before = std::max(cv.press_before, cv.press_floor);
//   cv.press = cv.press_before;

//   // Set temp_mu and temp
//   cv.temp_mu = cv.press/cv.den*punit->Temperature_mu;
//   cv.mu = muH_/(1.1 + cv.x_e - cv.x_h2);
//   cv.temp  = cv.mu*cv.temp_mu;
//   cv.temp_ = cv.temp*(1.0 + dlntemp_);

//   // SetRadiationVariables(pmb, cv);
//   cv.chi_pe = chi_pe;
//   cv.chi_ci = chi_ci;
//   cv.chi_co = chi_co;
  
//   cv.xi_ph_hi = xi_ph_hi;
//   cv.xi_ph_h2 = xi_ph_h2;
//   cv.xi_diss_h2 = xi_diss_h2;

//   // SetCosmicRayVariables(pmb, cv);
//   cv.xi_cr = xi_cr;

//   return;
// }

void CoolingSolverTigress::CoolingExplicitSubcycling(CoolVar& cv, const Real t_end) {

  int flag_bad_dt_sub;
  int nbad_dt_max = 3;
  int nsub = 0;
  Real t_done = 0.0, t_left = t_end, dt_sub;
  
  while ((nsub < nsub_max_) && (t_done < t_end)) {
    dt_sub = std::max(0.0, t_left);
    flag_bad_dt_sub = DoOneSubstep(cv, dt_sub);
    nsub++;
    t_done += dt_sub;
    if (flag_bad_dt_sub > nbad_dt_max) {
      // If a cell gets flagged repeately, it's usually a bad cell produced by
      // Riemann solver and has extremely low pressure. There is not much
      // cooling solver can do to fix this.
      t_done = t_end;
      // output_log_too_many_bad_dt(flag_bad_dt_sub);
    }
  }
  
  std::cout << "nH " << cv.nH << " nsub " << nsub
            << " nbad_dt " << flag_bad_dt_sub << std::endl;
  
  return;
}

int CoolingSolverTigress::DoOneSubstep(CoolVar& cv, Real& dt_sub) {

  // Following the practice adopated Anninos et al. (1997) and others,
  // sequentially update temperature followed by species abundances (H2 and HII
  // only). Alternatively, one can update chemistry first, although it's likely
  // Anninos et al. tried both and adopted the temperature-first method for a
  // reason.

  const int nbad_dt_max = 3;
  int flag_bad_dt_sub = 0;

  // Update dt_sub if 
  // cfl_cool_sub * MIN(|t_cool_net|, |t_chem|)
  // is smaller than the current estimate.
  // cfl_cool_sub is usually set to 0.1 (10% rule)
  CalculateCoolingRates(cv);
  CalculateChemicalRates(cv);
  dt_sub = std::min(dt_sub, cfl_cool_sub_/std::fabs(cv.it_cool_net));
  dt_sub = std::min(dt_sub, cfl_cool_sub_/cv.it_hii);
  if (cv.it_h2 != HUGE_NUMBER)
    // Skipped for hot gas
    dt_sub = std::min(dt_sub, cfl_cool_sub_/cv.it_h2);
  
  while (1) {
    int temp_ok = UpdateTemperature(cv, dt_sub);
    if (temp_ok) {
      // Recompute rate coefficients with updated T
      CalculateChemicalRates(cv);
      UpdateChemistry(cv, dt_sub);
      if (cv.len_shld > 0.0) {
        UpdateRadiationField(cv);
      }
      break;
    } else {
      // Real temp_mu_next = cv.temp_mu + cv.dtemp_mu_heat - cv.dtemp_mu_cool;
      // if ((temp_mu_next < 0) && (flag_bad_dt_sub > 1)) {
        // output_log_negative_temp_mu
      // }
      dt_sub *= 0.5;
      flag_bad_dt_sub++;
      if (flag_bad_dt_sub > nbad_dt_max) break;
    }
  }
  return flag_bad_dt_sub;
}

void CoolingSolverTigress::CalculateCoolingRates(CoolVar& cv) {

  // minimum heating rate coefficient [erg/s/H]
  const Real heat_min = 1.0e-35;
  
  Real cool_hyd = 0.0, cool_hyd_ = 0.0;
  Real cool_other = 0.0, cool_other_ = 0.0;
  Real cool_hyd_tbl = 0.0, cool_hyd_tbl_ = 0.0;
  Real cool_he_tbl = 0.0, cool_he_tbl_ = 0.0;
  Real cool_metal_tbl = 0.0, cool_metal_tbl_ = 0.0;
  Real cool = 0.0, cool_= 0.0;
  Real heat = 0.0, heat_= 0.0;
  
  // // Variables related to dust cooling
  // Real Tdust_in = pCVar->Tdust_in;
  // Real heating_dustUV = pCVar->heating_dustUV;
  Real cool_dust_net = 0.0, cool_dust_net_ = 0.0;
  Real temp_dust_eq, temp_dust_eq_;

  // Hydrogen cooling (assume that abundances remain fixed)
  // apply H species abundance cuts

  if (cv.x_hi > x_hi_cut) {
    cool_hyd  += CoolingHI(cv.x_e, cv.x_hi, cv.nH, cv.temp );
    cool_hyd_ += CoolingHI(cv.x_e, cv.x_hi, cv.nH, cv.temp_);
  }
  
  if (cv.x_hii > x_hii_cut) {
    cool_hyd  += CoolingHII(cv.x_e, cv.x_hii, cv.nH, cv.temp );
    cool_hyd_ += CoolingHII(cv.x_e, cv.x_hii, cv.nH, cv.temp_);
  }
  
  if ((cv.x_h2 > x_h2_cut) && (cv.temp < temp_hot1_)) {
    if (flag_h2_cooling) {
      cool_hyd  += CoolingH2(cv.x_h2, cv.x_hi, cv.nH, cv.temp );
      cool_hyd_ += CoolingH2(cv.x_h2, cv.x_hi, cv.nH, cv.temp_);
    }
    if (flag_h2_heating) {
    heat  += HeatingH2(cv.x_hi, cv.x_h2, cv.nH, cv.temp ,
                       cv.xi_diss_h2, cv.z_d);
    heat_ += HeatingH2(cv.x_hi, cv.x_h2, cv.nH, cv.temp_,
                       cv.xi_diss_h2, cv.z_d);
    }
  }

  // Other heating and cooling
  if (cv.temp < temp_hot1_) {
    cool_other  = CoolingOther(cv.x_e, cv.x_hi, cv.x_h2, cv.x_hii,
                               cv.nH, cv.temp , cv.dvdr, cv.z_d, cv.z_g,
                               cv.xi_cr, cv.chi_fuv, cv.chi_ci, cv.chi_co);
    cool_other_ = CoolingOther(cv.x_e, cv.x_hi, cv.x_h2, cv.x_hii,
                               cv.nH, cv.temp_, cv.dvdr, cv.z_d, cv.z_g,
                               cv.xi_cr, cv.chi_fuv, cv.chi_ci, cv.chi_co);
    heat  += HeatingOther(cv.x_e, cv.x_hi, cv.x_h2, cv.nH, cv.temp , cv.z_d, cv.xi_cr,
                          cv.chi_fuv, cv.xi_ph_h2, cv.xi_diss_h2, cv.xi_ph_hi);
    heat_ += HeatingOther(cv.x_e, cv.x_hi, cv.x_h2, cv.nH, cv.temp_, cv.z_d, cv.xi_cr,
                          cv.chi_fuv, cv.xi_ph_h2, cv.xi_diss_h2, cv.xi_ph_hi);
    
    // To save the cost, temp_dust_eq and associated cooling is computed only if temp <
    // temp_hot1. This implies hat temp_dust_eq in hot gas is from previous calculations
    // and is not used.
    
    if (flag_cool_dust_) {
      cool_dust_net  = get_cooling_dust_net(cv.temp_dust, cv.nH, cv.temp , cv.z_d,
                                            cv.heating_dust_uv, &temp_dust_eq );
      // Using the previous Td_eq as an initial guess may help the root-finding process
      cool_dust_net_ = get_cooling_dust_net(temp_dust_eq, cv.nH, cv.temp_, cv.z_d,
                                            cv.heating_dust_uv, &temp_dust_eq_);
      // Save equilibrium dust temperature to CVar
      cv.temp_dust = temp_dust_eq;
      if (cool_dust_net > 0.0) {
        cool_other  += cool_dust_net ;
        cool_other_ += cool_dust_net_;
      } else {
        heat  += cool_dust_net ;
        heat_ += cool_dust_net_;
      }
    }
  } else {
    heat  += heat_min;
    heat_ += heat_min;
  }

  if (cv.temp >= temp_hot0_) {
    cool_he_tbl  = cv.nH*Lambda_He_tbl(cv.temp );
    cool_he_tbl_ = cv.nH*Lambda_He_tbl(cv.temp_);
    cool_metal_tbl  = cv.nH*Lambda_metal_tbl(cv.temp , cv.z_g);
    cool_metal_tbl_ = cv.nH*Lambda_metal_tbl(cv.temp_, cv.z_g);
    if (flag_cool_hyd_cie_) {
      cool_hyd_tbl  = cv.nH*Lambda_H_tbl(cv.temp );
      cool_hyd_tbl_ = cv.nH*Lambda_H_tbl(cv.temp_);
    }
  }
  
  // Non-equilibrium hydrogen cooling
  cool  = cool_hyd ;
  cool_ = cool_hyd_;

  // Determine total cooling rate based on temperature
  if (cv.temp >= temp_hot0_) {
    if (cv.temp < temp_hot1_) {
      Real wgt2  = 1.0/(1.0 + std::exp(-10*(cv.temp  - 0.5*(temp_hot0_ + temp_hot1_))
                                       /(temp_hot1_ - temp_hot0_)));
      Real wgt2_ = 1.0/(1.0 + std::exp(-10*(cv.temp_ - 0.5*(temp_hot0_ + temp_hot1_))
                                       /(temp_hot1_ - temp_hot0_)));
      Real wgt1  = 1.0 - wgt2 ;
      Real wgt1_ = 1.0 - wgt2_;
      // smooth transition of hydrogen cooling from non-equilibrium to CIE
      if (flag_cool_hyd_cie_) {
        cool  = wgt1 *cool_hyd  +  wgt2 *cool_hyd_tbl ;
        cool_ = wgt1_*cool_hyd_ +  wgt2_*cool_hyd_tbl_;
      }
      cool  += wgt1 *cool_other  + cool_he_tbl  + wgt2 *cool_metal_tbl ;
      cool_ += wgt1_*cool_other_ + cool_he_tbl_ + wgt2_*cool_metal_tbl_;
      heat  *= wgt1 ;
      heat_ *= wgt1_;
    } else {
      // overwrite hydrogen cooling with CIE
      if (flag_cool_hyd_cie_) {
        cool  = cool_hyd_tbl ;
        cool_ = cool_hyd_tbl_;
      }
      cool  += cool_he_tbl  + cool_metal_tbl ;
      cool_ += cool_he_tbl_ + cool_metal_tbl_;
    }
  } else {
    cool  += cool_other ;
    cool_ += cool_other_;
  }
  
  // Cooling and heating rates per unit volume in code unit
  cv.cool_rate = (cv.nH*cool)*(punit->Time/punit->EnergyDensity);
  cv.heat_rate = (cv.nH*heat)*(punit->Time/punit->EnergyDensity);

  // Derivative of net cooling rate in code unit
  cv.df_dtemp_mu = cv.mu*cv.nH*((heat_ - heat) -
                                (cool_ - cool))/(cv.temp_ - cv.temp);
  cv.df_dtemp_mu *= (punit->Time/punit->EnergyDensity)*punit->Temperature_mu;
  
  // inverse of cooling/heating time in code unit
  cv.it_cool = cv.cool_rate/(igm1_*cv.press);
  cv.it_heat = cv.heat_rate/(igm1_*cv.press);
  cv.it_cool_net = cv.it_cool - cv.it_heat;

  return;

}

void CoolingSolverTigress::CalculateChemicalRates(CoolVar& cv) {

  // ionized hydrogen
  HIIRates(cv.nH, cv.temp, cv.x_hi, cv.x_h2, cv.x_e,
           cv.z_d, cv.xi_ph_hi, cv.xi_cr, cv.chi_fuv,
           &(cv.crate_hii), &(cv.drate_hii));
  cv.crate_hii *= punit->Time;
  cv.drate_hii *= punit->Time;
  cv.it_hii = std::fabs(cv.x_hi*cv.crate_hii - cv.x_hii*cv.drate_hii);

  // molecular hydrogen  
  if (cv.temp >= temp_hot1_) {
    // Do not consider chemical timescale in substepping
    cv.crate_h2 = 0.0;
    cv.drate_h2 = 0.0;
    cv.it_h2 = HUGE_NUMBER;
  } else {
    H2Rates(cv.nH, cv.temp, cv.x_hi, cv.x_h2, cv.x_e, cv.z_d,
            cv.xi_ph_h2, cv.xi_diss_h2, cv.xi_cr,
            &(cv.crate_h2), &(cv.drate_h2));
    cv.crate_h2 *= punit->Time;
    cv.drate_h2 *= punit->Time;
    cv.it_h2 = std::fabs(cv.x_hi*cv.crate_h2 - cv.x_h2*cv.drate_h2);
  }
  
  return;
}

int CoolingSolverTigress::UpdateTemperature(CoolVar& cv, const Real dt) {
  
  // In Eq. 59 in Kim, J.-G. et al. (2023) The term involving the derivative
  // reduces to -(gamma-1)*mH*(Delta t)/(rho*kB)*(dH/dT_mu)
  Real deriv_term = -cv.df_dtemp_mu/(igm1_*cv.den)*dt;
  Real dtemp_mu_heat_noderiv = cv.temp_mu*cv.it_heat*dt;
  Real dtemp_mu_cool_noderiv = cv.temp_mu*cv.it_cool*dt;
  cv.dtemp_mu_heat = dtemp_mu_heat_noderiv/(1.0 + deriv_term);
  cv.dtemp_mu_cool = dtemp_mu_cool_noderiv/(1.0 + deriv_term);

  Real temp_mu_next = cv.temp_mu + cv.dtemp_mu_heat - cv.dtemp_mu_cool;

  if ((temp_mu_next < 0) || (std::isnan(temp_mu_next)) ||
      (temp_mu_next > (1.0 + 2.0*cfl_cool_sub_)*cv.temp_mu) ||
      (temp_mu_next < (1.0 - 2.0*cfl_cool_sub_)*cv.temp_mu)) {
    // Leave message for debugging purpose
    //
    return 0;
  } else {
    // Successful update
    Real press_next = cv.den*temp_mu_next/punit->Temperature_mu;
    // heating due to flooring
    if (press_next < cv.press_floor) {
      cv.dpress_floor += cv.press_floor - press_next;
      cv.press = cv.press_floor;
    } else {
      cv.press = press_next;
    }
  }
  
  // Accumulate pressure change by cooling and heating separately
  cv.dpress_cool += cv.den*cv.dtemp_mu_cool/punit->Temperature_mu;
  cv.dpress_heat += cv.den*cv.dtemp_mu_heat/punit->Temperature_mu;

  // Set temp_mu and temp using updated values
  cv.temp_mu = cv.press/cv.den*punit->Temperature_mu;
  cv.mu = muH_/(1.1 + cv.x_e - cv.x_h2);
  cv.temp  = cv.mu*cv.temp_mu;
  cv.temp_ = cv.temp*(1.0 + dlntemp_);
  
  return 1;
}

void CoolingSolverTigress::UpdateChemistry(CoolVar& cv, const Real dt) {

  // Update species abundances semi-implicitly
  // We write the ODE for xH2 as
  // d xH2/dt  = C_H2 *xHI - D_H2 *xH2
  // d xHII/dt = C_HII*xHI - D_HII*xH2
  // where xHI = 1.0 - 2*xH2 - xHII
  // We (incorrectly) assume that the above equations are linear
  // in xH2 and xHII.
  // Using the backward differentiation formula
  // we can explictly solve for xH2^(n+1) and xHII^(n+1)

  Real x_hii, x_h2, x_hi, x_e, x_cii;

  Real a = 1.0 + (2.0*cv.crate_h2 + cv.drate_h2)*dt;
  Real b = cv.crate_h2*dt;
  Real c = 2.0*cv.crate_hii*dt;
  Real d = 1.0 + (cv.crate_hii + cv.drate_hii)*dt;
  Real e = cv.x_h2 + cv.crate_h2*dt;
  Real f = cv.x_hii + cv.crate_hii*dt;

  x_hii = (a*f - c*e)/(a*d - b*c);
  x_hii = std::min(1.0, x_hii);
  x_hii = std::max(TINY_NUMBER, x_hii);

  if (cv.temp < temp_hot1_) {
    x_h2 = (d*e - b*f)/(a*d - b*c);
  } else {
    x_h2 = 0.0; // no H2 in hot gas
  }

  // Enforce closure relation (similar to approach taken by RAMSES)
  Real x_sum = x_hii + 2.0*x_h2;
  if (x_sum > 1.0) {
    x_hii = x_hii/x_sum;
    x_h2 = x_h2/x_sum;
  }
  if (x_h2 < 0.25) {
    x_h2 = std::max(TINY_NUMBER, x_h2);
    x_h2 = std::min(0.5, x_h2);
    x_hii = std::max(TINY_NUMBER, x_hii);
    x_hii = std::min(1.0 - 2.0*x_h2, x_hii);
  } else {
    x_hii = std::max(TINY_NUMBER, x_hii);
    x_hii = std::min(1.0, x_hii);
    x_h2 = std::max(TINY_NUMBER, x_h2);
    x_h2 = std::min(1.0 - x_hii, x_h2);
  }
  x_hi = std::max(1.0 - x_hii - 2.0*x_h2, TINY_NUMBER);

  cv.x_h2 = x_h2;
  cv.x_hii = x_hii;
  cv.x_hi = x_hi;

  cv.x_e = CalculateFreeElectron(cv.x_e, cv.x_hi, cv.x_h2, cv.x_hii,
                                 cv.nH, cv.temp, cv.z_d, cv.z_g, cv.xi_cr,
                                 cv.chi_fuv, cv.chi_ci, &x_cii);
  
  Real x_oii = cv.x_hii*xOstd*cv.z_g;
  Real x_co = CalculateCOAbundance(cv.x_h2, x_cii, x_oii, cv.nH, cv.z_d, cv.z_g,
                                   cv.xi_cr, cv.chi_co);
  cv.x_ci = std::max(xCstd*cv.z_g - x_cii - x_co, TINY_NUMBER);
  
  return;
}

void CoolingSolverTigress::UpdateRadiationField(CoolVar& cv) {

  
  // Update H2-dissociating, CI-ionizing radiation fields based on updated abundances
  Real col_shld_h2_cgs = 0.5*(cv.nH*cv.x_h2)*(cv.len_shld*punit->Length);

  // Calculate shielding factors
  Real x = col_shld_h2_cgs/5.0e14;
  Real sqrt_onepx = std::sqrt(1 + x);
  Real fshld_h2 = 0.965/SQR((1.0 + x*b5_inv)) +
    0.034970262640168365/sqrt_onepx*exp(-8.5e-4*sqrt_onepx);

  Real col_shld_ci_cgs = 0.5*(cv.nH*cv.x_ci)*(cv.len_shld*punit->Length);
  Real tau_ci = 1.6e-17*col_shld_ci_cgs;
  Real r_h2 = x*1.4e-7;
  Real fshld_ci = exp(-tau_ci)*exp(-r_h2)/(1.0 + r_h2);
  
  cv.chi_h2 = cv.chi_lw*fshld_h2;
  cv.xi_diss_h2 = cv.chi_h2*xi_diss_h2_isrf;
  cv.chi_ci = cv.chi_lw*fshld_ci;
  
}

// void CoolingSolverTigress::UpdateHydroVariables(MeshBlock *pmb, CoolVar& cv,
//                                                 const int k, const int j, const int i) {

//   //Real dt_mhd = pmb->pmy_mesh->dt;
//   Real dt_mhd = 1.0;
//   Real gamma = 5.0/3.0;
//   Real igm1_ = 1.0/(gamma - 1.0);
//   Real igm1_idt_mhd = 1/(gamma - 1.0)/dt_mhd;

//   // save instantaneous rates at the end of subcycling;
//   // cv.dpress_heat = cv.press*cv.it_heat*t_end;
//   // cv.dpress_cool = cv.press*cv.it_cool*t_end;

//   // Save pressure difference due to cooling
//   // Real dpress = cv.press - cv.press_before;
  
//   // Apply presure floor after solving cooling
//   if (cv.press < cv.press_floor)
//     cv.dpress_floor += cv.press_floor - cv.press;
//   cv.press = std::max(cv.press, cv.press_floor);

//   // if (bookkeeping) {
//   //   edot(k,j,i) = dpress*igm1_idt_mhd;
//   //   edot_floor(k,j,i) = delta_press_floor*igm1_idt_mhd;
//   // }

//   // change internal energy
//   // pmb->phydro->u(IEN,k,j,i) = cv.press*igm1_ + cv.e_non_thermal;
//   // pmb->phydro->w(IPR,k,j,i) = cv.press;

//   return;

// }

void CoolingSolverTigress::UpdateVariables1D(Real *den, Real *press,
                                             Real *x_h2, Real *x_hii, Real *x_e,
                                             Real *chi_h2, Real *chi_ci,
                                             CoolVar& cv, const int i) {
  
  //Real dt_mhd = pmb->pmy_mesh->dt;
  //  Real dt_mhd = 1.0;
  
  // Real igm1_ = 1.0/(gamma - 1.0);
  // Real igm1_idt_mhd = 1/(gamma - 1.0)/dt_mhd;

  // save instantaneous rates at the end of subcycling;
  // cv.dpress_heat = cv.press*cv.it_heat*t_end;
  // cv.dpress_cool = cv.press*cv.it_cool*t_end;

  // Save pressure difference due to cooling
  // Real dpress = cv.press - cv.press_before;
  
  // Apply presure floor after solving cooling
  if (cv.press < cv.press_floor)
    cv.dpress_floor += cv.press_floor - cv.press;
  cv.press = std::max(cv.press, cv.press_floor);

  // if (bookkeeping) {
  //   edot(k,j,i) = dpress*igm1_idt_mhd;
  //   edot_floor(k,j,i) = delta_press_floor*igm1_idt_mhd;
  // }

  // change internal energy
  // den[i] = cv.press*igm1_ + cv.e_non_thermal;

  press[i] = cv.press;
  x_h2[i] = cv.x_h2;
  x_hii[i] = cv.x_hii;
  x_e[i] = cv.x_e;
  chi_h2[i] = cv.chi_h2;
  chi_ci[i] = cv.chi_ci;
  
  return;

}
