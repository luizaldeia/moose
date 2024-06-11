//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TwoPhaseFluidProperties.h"
#include "SinglePhaseFluidProperties.h"
#include "Numerics.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

// class TwoPhaseFluidProperties;

/**
 * Class for the Homogeneous Equilibrium Model fluid properties as a function of pressure and
 * temperature for water
 *
 *
 */
class HEM : public TwoPhaseFluidProperties
{
protected:
  ///  Loading the 2-phase fluid properties UserObject
  const TwoPhaseFluidProperties & _fp_2phase;
  /// pointer to liquid fluid properties object (if provided 2-phase object)
  const SinglePhaseFluidProperties & _fp_liquid;
  /// pointer to vapor fluid properties object (if provided 2-phase object)
  const SinglePhaseFluidProperties & _fp_vapor;
  /// Initial temperature guess for the fluid state classification
  const Real & _temp_guess;

public:
  static InputParameters validParams();

  HEM(const InputParameters & parameters);

  /**
   * Returns the critical pressure
   */
  virtual Real p_critical() const override { return _fp_2phase.p_critical(); }

  /**
   * Computes the saturation temperature at a pressure
   *
   * @param[in] p  pressure
   */
  virtual Real T_sat(Real p) const override { return _fp_2phase.T_sat(p); }
  virtual ADReal T_sat(const ADReal & p) const override { return _fp_2phase.T_sat(p); }

  /**
   * Computes the saturation pressure at a temperature
   *
   * @param[in] T  temperature
   */
  virtual Real p_sat(Real T) const override { return _fp_2phase.p_sat(T); }
  virtual ADReal p_sat(const ADReal & T) const override { return _fp_2phase.p_sat(T); }

  /**
   * Computes dT/dp along the saturation line
   *
   * @param[in] p  pressure
   */
  virtual Real dT_sat_dp(Real p) const override { return _fp_2phase.dT_sat_dp(p); }
  // virtual ADReal dT_sat_dp(const ADReal & p) const override { return _fp_2phase.dT_sat_dp(p); }

  /**
   * Computes latent heat of vaporization
   *
   * @param p  pressure
   * @param T  temperature
   */
  virtual Real h_lat(Real p, Real T) const override { return _fp_2phase.h_lat(p, T); }
  virtual ADReal h_lat(const ADReal & p, const ADReal & T) const override
  {
    return _fp_2phase.h_lat(p, T);
  }

  // virtual bool supportsPhaseChange() const override { return true; }
  virtual bool supportsPhaseChange() const override { return _fp_2phase.supportsPhaseChange(); }

  /**
   * Mixture density from pressure and temperature
   *
   * @param[in] p       pressure
   * @param[in] T       temperature
   * @param[in] alpha   void fraction
   * @return            mixture density
   */
  virtual ADReal rho_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const;

  /**
   * Liquid phase density from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        liquid phase density
   */
  virtual ADReal rho_liquid_from_p_T(ADReal p, ADReal T) const;

  /**
   * Vapor phase density from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        vapor phase density
   */
  virtual ADReal rho_vapor_from_p_T(ADReal p, ADReal T) const;

  /**
   * Mixture specific volume from pressure and temperature
   *
   * @param[in] p       pressure
   * @param[in] T       temperature
   * @param[in] alpha   void fraction
   * @return            mixture density
   */
  virtual ADReal v_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const;

  /**
   * Liquid phase specific volume from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        liquid phase density
   */
  virtual ADReal v_liquid_from_p_T(ADReal p, ADReal T) const;

  /**
   * Vapor phase specific volume from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        liquid phase density
   */
  virtual ADReal v_vapor_from_p_T(ADReal p, ADReal T) const;

  /**
   * Mixture internal energy from pressure and temperature
   *
   * @param[in] p       pressure
   * @param[in] T       temperature
   * @param[in] alpha   void fraction
   * @return            mixture density
   */
  virtual ADReal e_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const;

  /**
   * Liquid phase internal energy from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        liquid phase density
   */
  virtual ADReal e_liquid_from_p_T(ADReal p, ADReal T) const;

  /**
   * Vapor phase internal energy from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        vapor phase density
   */
  virtual ADReal e_vapor_from_p_T(ADReal p, ADReal T) const;

  /**
   * Liquid phase internal energy from pressure and density
   *
   * @param[in] p     pressure
   * @param[in] rho   density
   * @return          liquid phase density
   */
  virtual ADReal e_liquid_from_p_rho(ADReal p, ADReal rho) const;

  /**
   * Vapor phase internal energy from pressure and density
   *
   * @param[in] p     pressure
   * @param[in] rho   density
   * @return          vapor phase density
   */
  virtual ADReal e_vapor_from_p_rho(ADReal p, ADReal rho) const;

  /**
   * Mixture internal energy from pressure and density
   *
   * @param[in] p       pressure
   * @param[in] rho     density
   * @param[in] alpha   void fraction
   * @return            Mixture density
   */
  virtual ADReal e_mixture_from_p_rho(ADReal p, ADReal rho, ADReal alpha) const;

  /**
   * Mixture specific enthalpy from pressure and temperature
   *
   * @param[in] p       pressure
   * @param[in] T       temperature
   * @param[in] alpha   void fraction
   * @return            mixture density
   */
  virtual ADReal h_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const;

  /**
   * Liquid phase specific enthalpy from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        liquid phase density
   */
  virtual ADReal h_liquid_from_p_T(ADReal p, ADReal T) const;

  /**
   * Vapor phase specific enthalpy from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        vapor phase density
   */
  virtual ADReal h_vapor_from_p_T(ADReal p, ADReal T) const;

  /**
   * Mixture specific entropy from pressure and temperature
   *
   * @param[in] p       pressure
   * @param[in] T       temperature
   * @param[in] alpha   void fraction
   * @return            mixture density
   */
  virtual ADReal s_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const;

  /**
   * Liquid phase specific entropy from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        liquid phase density
   */
  virtual ADReal s_liquid_from_p_T(ADReal p, ADReal T) const;

  /**
   * Vapor phase specific entropy from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        vapor phase density
   */
  virtual ADReal s_vapor_from_p_T(ADReal p, ADReal T) const;

  /**
   * Mixture thermal conductivity from pressure and temperature
   *
   * @param[in] p       pressure
   * @param[in] T       temperature
   * @param[in] alpha   void fraction
   * @return            mixture thermal conductivity
   */
  virtual ADReal k_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const;

  /**
   * Liquid phase thermal conductivity from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        liquid phase thermal conductivity
   */
  virtual ADReal k_liquid_from_p_T(ADReal p, ADReal T) const;

  /**
   * Vapor phase thermal conductivity from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        vapor phase thermal conductivity
   */
  virtual ADReal k_vapor_from_p_T(ADReal p, ADReal T) const;

  /**
   * Mixture constant-pressure specific heat from pressure and temperature
   *
   * @param[in] p       pressure
   * @param[in] T       temperature
   * @param[in] alpha   void fraction
   * @return            mixture constant-pressure specific heat
   */
  virtual ADReal cp_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const;

  /**
   * Liquid phase constant-pressure specific heat from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        liquid phase constant-pressure specific heat
   */
  virtual ADReal cp_liquid_from_p_T(ADReal p, ADReal T) const;

  /**
   * Vapor phase constant-pressure specific heat from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        vapor phase constant-pressure specific heat
   */
  virtual ADReal cp_vapor_from_p_T(ADReal p, ADReal T) const;

  /**
   * Mixture constant-volume specific heat from pressure and temperature
   *
   * @param[in] p       pressure
   * @param[in] T       temperature
   * @param[in] alpha   void fraction
   * @return            mixture constant-pressure specific heat
   */
  virtual ADReal cv_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const;

  /**
   * Liquid phase constant-volume specific heat from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        liquid phase constant-pressure specific heat
   */
  virtual ADReal cv_liquid_from_p_T(ADReal p, ADReal T) const;

  /**
   * Vapor phase constant-volume specific heat from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        vapor phase constant-pressure specific heat
   */
  virtual ADReal cv_vapor_from_p_T(ADReal p, ADReal T) const;

  /**
   * Mixture speed of sound from pressure and temperature
   *
   * @param[in] p       pressure
   * @param[in] T       temperature
   * @param[in] alpha   void fraction
   * @return        mixture speed of sound
   */
  virtual ADReal c_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const;

  /**
   * Liquid phase speed of sound from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        liquid phase speed of sound
   */
  virtual ADReal c_liquid_from_p_T(ADReal p, ADReal T) const;

  /**
   * Vapor phase speed of sound from pressure and temperature
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        vapor phase speed of sound
   */
  virtual ADReal c_vapor_from_p_T(ADReal p, ADReal T) const;

  /**
   * Liquid phase speed of sound from specific volume and internal energy
   *
   * @param[in] v   specific volume
   * @param[in] e   internal energy
   * @return        liquid phase speed of sound
   */
  virtual ADReal c_liquid_from_v_e(ADReal v, ADReal e) const;

  /**
   * Vapor phase speed of sound from specific volume and internal energy
   *
   * @param[in] v   specific volume
   * @param[in] e   internal energy
   * @return        vapor phase speed of sound
   */
  virtual ADReal c_vapor_from_v_e(ADReal v, ADReal e) const;

  /**
   * Liquid phase pressure from specific volume and internal energy
   *
   * @param[in] v   specific volume
   * @param[in] e   internal energy
   * @return        liquid phase pressure
   */
  virtual ADReal p_liquid_from_v_e(ADReal v, ADReal e) const;

  /**
   * Vapor phase pressure from specific volume and internal energy
   *
   * @param[in] v   specific volume
   * @param[in] e   internal energy
   * @return        vapor phase pressure
   */
  virtual ADReal p_vapor_from_v_e(ADReal v, ADReal e) const;

  /**
   * Liquid phase temperature from specific volume and internal energy
   *
   * @param[in] v   specific volume
   * @param[in] e   internal energy
   * @return        liquid phase temperature
   */
  virtual ADReal T_liquid_from_v_e(ADReal v, ADReal e) const;

  /**
   * Vapor phase temperature from specific volume and internal energy
   *
   * @param[in] v   specific volume
   * @param[in] e   internal energy
   * @return        vapor phase temperature
   */
  virtual ADReal T_vapor_from_v_e(ADReal v, ADReal e) const;

  /**
   * Mixture dynamic viscosity from specific volume and internal energy
   *
   * @param[in] p   pressure
   * @param[in] T   temperature
   * @return        Mixture dynamic viscosity
   */
  virtual ADReal mu_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const;

  /**
   * Liquid phase dynamic viscosity from specific volume and internal energy
   *
   * @param[in] v   specific volume
   * @param[in] e   internal energy
   * @return        liquid phase dynamic viscosity
   */
  virtual ADReal mu_liquid_from_v_e(ADReal v, ADReal e) const;

  /**
   * Vapor phase dynamic viscosity from specific volume and internal energy
   *
   * @param[in] v   specific volume
   * @param[in] e   internal energy
   * @return        vapor phase dynamic viscosity
   */
  virtual ADReal mu_vapor_from_v_e(ADReal v, ADReal e) const;

  // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  // Calculating the relevant fluid properties as a function of the conservative variables:
  // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

  /**
   * Fluid internal energy from rho*E*A, rho*u*A, rho*A, and A: rhoe = rhoEA/A -
   * 0.5*rho*vel^2
   *
   * @param[in] rhoEA  energy conservative variable
   * @param[in] rhouA  momentum conservative variable
   * @param[in] rhoA   mass conservative variable
   * @param[in] A      cross-sectional area
   * @return           Fluid density times internal energy
   */
  virtual ADReal
  rhoe_from_rhoEA_rhouA_rhoA_A(ADReal rhoEA, ADReal rhouA, ADReal rhoA, ADReal A) const;

  // determine the fluid state from rho*E*A, rho*u*A, rho*A, and A
  struct HEMState
  {
    ADReal rho;
    ADReal v;
    ADReal vel;
    ADReal e;
    ADReal p;
    ADReal T;
    ADReal h;
    ADReal H;
    ADReal c;
    ADReal cp;
    ADReal cv;
    ADReal k;
    ADReal alpha;
  };
  /**
   * Fluid state
   *
   * @param[in] rhoEA           energy conservative variable
   * @param[in] rhouA           momentum conservative variable
   * @param[in] rhoA            mass conservative variable
   * @param[in] A               cross-sectional area
   * @param[in] temp_guess      Temperature initial guess
   * @return                    Fluid state
   */
  virtual HEMState fluid_state(ADReal rhoA, ADReal rhouA, ADReal rhoEA, ADReal A) const;
};

#pragma GCC diagnostic pop
