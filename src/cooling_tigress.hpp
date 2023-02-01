#ifndef MICROPHYSICS_COOLING_TIGRESS_HPP_
#define MICROPHYSICS_COOLING_TIGRESS_HPP_

#include "defs.hpp"
#include "units.hpp"
#include "parameter_input.hpp"

// Cooling-related POD data holder
struct CoolVar {

  Real z_g, z_d;
  Real den, nH;
  Real press, press_before, press_floor;
  Real e_non_thermal;

  // Note that 
  // P = (gamma-1) e_thermal = rho * kB * temp_mu / mH
  // where temp_mu = temp / mu
  // and mu = muH / (1.1 + x_e - x_H2)
  Real temp_mu, mu;
  Real temp, temp_;

  Real dtemp_mu_cool, dtemp_mu_heat;
  Real dpress_heat, dpress_cool, dpress_floor;

  // Abundances
  // Define x_h2 = n(H2)/nH so that 0 < x_h2 < 0.5
  Real x_e, x_h2, x_hi, x_hii, x_ci;

  // local radiation field and cosmic-ray ionization rate
  Real chi_fuv, chi_lw, chi_h2, chi_ci, chi_co;
  Real xi_ph_hi, xi_ph_h2, xi_diss_h2;
  Real xi_cr;

  Real dvdr;

  // shielding variables
  Real len_shld;
  
  // cooling and heating rates, derivate of net cooling, and inverse of
  // timescales
  Real cool_rate, heat_rate;
  Real df_dtemp_mu;
  Real it_cool, it_heat, it_cool_net;
  
  // chemical rates and inverse of timescales
  Real crate_hii, drate_hii;
  Real crate_h2, drate_h2;
  Real it_hii, it_h2;

  Real heating_dust_uv;
  Real temp_dust;
  
};

class CoolingSolverTigress {
 public:
  explicit CoolingSolverTigress(ParameterInput *pin, Units *punit);
  ~CoolingSolverTigress();
  
  //void OperatorSplitSolver(Real dt);
  
  void OperatorSplitSolver1D(Real *den, Real *press,
                             Real *x_h2, Real *x_hii, Real *x_e,
                             Real *chi_pe, Real *chi_lw,
                             Real *chi_h2, Real *chi_ci,
                             Real *xi_cr, Real *len_shld,
                             const Real chi0,
                             const Real z_g, const Real z_d,
                             int il, int iu, const Real dt);
  
  Units *punit;
  
 private:
  // void SetupCoolingVariables(MeshBlock *pmb, CoolVar& cv,
  //                            const int k, const int j, const int i);
  void SetupVariables1D(Real *den, Real *press,
                        Real *x_h2, Real *x_hii, Real *x_e,
                        Real *chi_pe, Real *chi_lw,
                        Real *chi_h2, Real *chi_ci,
                        Real *xi_cr, Real *len_shld,
                        const Real chi0,
                        const Real z_g, const Real z_d,
                        CoolVar& cv,
                        const int i);

  // void UpdateHydroVariables(MeshBlock *pmb, CoolVar& cv,
  //                           const int k, const int j, const int i);
  void UpdateVariables1D(Real *den, Real *press,
                         Real *x_h2, Real *x_hii, Real *x_e,
                         Real *chi_h2, Real *chi_ci,
                         CoolVar& cv, const int i);
  
  void CoolingExplicitSubcycling(CoolVar& cv, const Real t_end);
  int DoOneSubstep(CoolVar& cv, Real& dt_sub);
  void CalculateCoolingRates(CoolVar& cv);
  void CalculateChemicalRates(CoolVar& cv);
  int UpdateTemperature(CoolVar& cv, const Real dt);
  void UpdateChemistry(CoolVar& cv, const Real dt);
  void UpdateRadiationField(CoolVar& cv);

  const int flag_cool_dust_;
  const int flag_cool_hyd_cie_;
  const Real cfl_cool_sub_;
  const Real muH_;
  const Real sigma_dust_pe0_, sigma_dust_lw0_;
  const int nsub_max_;
  const Real code_den_to_nH_;
  const Real temp_hot0_, temp_hot1_;
  
  // relative change in temperature (used for derivative calculation)
  const Real dlntemp_ = 0.02;
  const Real gm1_, igm1_;
};

void HIIRates(const Real nH, const Real T, const Real xHI,
              const Real xH2, const Real xe, const Real Z_d,
              const Real xi_ph_HI, const Real xi_CR, const Real chi_fuv,
              Real *crate_HII, Real *drate_HII);

void H2Rates(const Real nH, const Real T, const Real xHI,
             const Real xH2, const Real xe, const Real Z_d,
             const Real xi_ph_H2, const Real xi_diss_H2, const Real xi_CR,
             Real *crate_H2, Real *drate_H2);
  
Real HeatingH2(const Real xHI, const Real xH2, const Real nH,
               const Real T, const Real xi_diss_H2, const Real Z_d);

Real HeatingOther(const Real xe, const Real xHI, const Real xH2,
                  const Real nH, const Real T, const Real Z_d,
                  const Real xi_CR, const Real chi_fuv,
                  const Real xi_ph_H2, const Real xi_diss_H2, const Real xi_ph_HI);

Real CoolingH2(const Real xH2, const Real xHI, const Real nH, const Real T);

Real CoolingHI(const Real xe, const Real xHI, const Real nH, const Real T);

Real CoolingHII(const Real xe, const Real xHII, const Real nH, const Real T);

Real CoolingOther(const Real xe, const Real xHI, const Real xH2, const Real xHII,
                  const Real nH, const Real T, const Real dvdr,
                  const Real Z_d, const Real Z_g, const Real xi_CR,
                  const Real chi_fuv, const Real chi_ci, const Real chi_co);

// Cooling function from CIE table
Real Lambda_H_tbl(const Real T);

Real Lambda_He_tbl(const Real T);

Real Lambda_metal_tbl(const Real T, const Real Z_g);

Real CalculateFreeElectron(const Real xe, const Real xHI, const Real xH2, const Real xHII,
                           const Real nH, const Real T, const Real Z_d, const Real Z_g,
                           const Real xi_CR, const Real chi_fuv, const Real chi_ci, Real *xCII_eq);



Real get_xH2(const Real nH, const Real T, const Real Z_d,
             const Real xi_CR, const Real chi_h2);

Real get_xe_fast(const Real nH, const Real T, const Real Z_d, const Real Z_g,
                 const Real xi_CR, const Real chi_fuv, const Real chi_ci, const Real chi_h2,
                 const Real xH2_in, const int iequil,
                 Real *pxHII, Real *pxH2, Real *pxCII, Real *pxCI);

Real CalculateCOAbundance(const Real xH2, const Real xCII, const Real xOII, const Real nH,
                          const Real Z_d, const Real Z_g, const Real xi_CR, const Real chi_co);

Real get_cooling_dust_net(const Real Td_in, const Real nH, const Real T, const Real Z_d,
                          const Real heatingdust_UV, Real *Td_out);

// Real CalculateSumofSqrt(Real a, Real b);
// Real HII(Real xe, Real xHII, Real nH, Real T);

#endif // MICROPHYSICS_COOLING_TIGRESS_HPP_
