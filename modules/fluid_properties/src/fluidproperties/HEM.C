//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HEM.h"

registerMooseObject("FluidPropertiesApp", HEM);

InputParameters
HEM::validParams()
{
  InputParameters params = TwoPhaseFluidProperties::validParams();
  params.addClassDescription("Class for the Homogeneous Equilibrium Model fluid properties as a "
                             "function of pressure and temperature for water");

  params.addRequiredParam<UserObjectName>("fp_2phase", "Two-phase water property user object name");

  // This is necessary because initialize() must be called before any interface
  // can be used (which can occur as early as initialization of variables).
  params.set<ExecFlagEnum>("execute_on") = EXEC_INITIAL;

  return params;
}

HEM::HEM(const InputParameters & parameters)
  : TwoPhaseFluidProperties(parameters),
    _fp_2phase(getUserObject<TwoPhaseFluidProperties>("fp_2phase")),
    _fp_liquid(getUserObjectByName<SinglePhaseFluidProperties>(_fp_2phase.getLiquidName())),
    _fp_vapor(getUserObjectByName<SinglePhaseFluidProperties>(_fp_2phase.getVaporName()))
{
}

// Definitions:
//  p      pressure [Pa]
//  T      temperature [K]
//  e      specific internal energy [J/kg]
//  v      specific volume [m^3/kg]
//  rho    density [kg/m^3]
//  h      specific enthalpy [J/kg]
//  s      specific entropy [J/(kg*K)]
//  mu     viscosity [Pa*s]
//  k      thermal conductivity [W/(m*K)]
//  c      speed of sound [m/s]
//  cp     constant-pressure specific heat [J/K]
//  cv     constant-volume specific heat [J/K]
//  beta   volumetric thermal expansion coefficient [1/K]
//  g      Gibbs free energy [J]
//  pp_sat partial pressure at saturation [Pa]
//  gamma  Adiabatic ratio (cp/cv) [-]

// Mixture density from pressure and temperature
ADReal
HEM::rho_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const
{
  return (1 - alpha) * (_fp_liquid.rho_from_p_T(p, T)) + alpha * (_fp_vapor.rho_from_p_T(p, T));
}

// Liquid phase density from pressure and temperature
ADReal
HEM::rho_liquid_from_p_T(ADReal p, ADReal T) const
{
  return _fp_liquid.rho_from_p_T(p, T);
}

// Vapor phase density from pressure and temperature
ADReal
HEM::rho_vapor_from_p_T(ADReal p, ADReal T) const
{
  return _fp_vapor.rho_from_p_T(p, T);
}

// Mixture specific volume from pressure and temperature
ADReal
HEM::v_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const
{
  // return (1 - alpha) / _fp_liquid.rho_from_p_T(p, T) + alpha / _fp_vapor.rho_from_p_T(p, T);
  return 1 / rho_mixture_from_p_T(p, T, alpha);
}

// Liquid phase specific volume from pressure and temperature
ADReal
HEM::v_liquid_from_p_T(ADReal p, ADReal T) const
{
  return 1 / _fp_liquid.rho_from_p_T(p, T);
}

// Vapor phase specific volume from pressure and temperature
ADReal
HEM::v_vapor_from_p_T(ADReal p, ADReal T) const
{
  return 1 / _fp_vapor.rho_from_p_T(p, T);
}

// Mixture internal energy from pressure and temperature
ADReal
HEM::e_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const
{
  return (1 - alpha) * (_fp_liquid.e_from_p_T(p, T)) + alpha * (_fp_vapor.e_from_p_T(p, T));
}

// Liquid phase internal energy from pressure and temperature
ADReal
HEM::e_liquid_from_p_T(ADReal p, ADReal T) const
{
  return _fp_liquid.e_from_p_T(p, T);
}

// Vapor phase internal energy from pressure and temperature
ADReal
HEM::e_vapor_from_p_T(ADReal p, ADReal T) const
{
  return _fp_vapor.e_from_p_T(p, T);
}

// Mixture specific enthalpy from pressure and temperature
ADReal
HEM::h_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const
{
  return (1 - alpha) * (_fp_liquid.h_from_p_T(p, T)) + alpha * (_fp_vapor.h_from_p_T(p, T));
}

// Liquid phase specific enthalpy from pressure and temperature
ADReal
HEM::h_liquid_from_p_T(ADReal p, ADReal T) const
{
  return _fp_liquid.h_from_p_T(p, T);
}

// Vapor phase specific enthalpy from pressure and temperature
ADReal
HEM::h_vapor_from_p_T(ADReal p, ADReal T) const
{
  return _fp_vapor.h_from_p_T(p, T);
}

// Mixture specific entropy from pressure and temperature
ADReal
HEM::s_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const
{
  return (1 - alpha) * (_fp_liquid.s_from_p_T(p, T)) + alpha * (_fp_vapor.s_from_p_T(p, T));
}

// Liquid phase specific entropy from pressure and temperature
ADReal
HEM::s_liquid_from_p_T(ADReal p, ADReal T) const
{
  return _fp_liquid.s_from_p_T(p, T);
}

// Vapor phase specific entropy from pressure and temperature
ADReal
HEM::s_vapor_from_p_T(ADReal p, ADReal T) const
{
  return _fp_vapor.s_from_p_T(p, T);
}

// Mixture thermal conductivity from pressure and temperature
ADReal
HEM::k_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const
{
  return (1 - alpha) * (_fp_liquid.k_from_p_T(p, T)) + alpha * (_fp_vapor.k_from_p_T(p, T));
}

// Liquid phase thermal conductivity from pressure and temperature
ADReal
HEM::k_liquid_from_p_T(ADReal p, ADReal T) const
{
  return _fp_liquid.k_from_p_T(p, T);
}

// Vapor phase thermal conductivity from pressure and temperature
ADReal
HEM::k_vapor_from_p_T(ADReal p, ADReal T) const
{
  return _fp_vapor.k_from_p_T(p, T);
}

// Mixture constant-pressure specific heat from pressure and temperature
ADReal
HEM::cp_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const
{
  return (1 - alpha) * (_fp_liquid.cp_from_p_T(p, T)) + alpha * (_fp_vapor.cp_from_p_T(p, T));
}

// Liquid phase constant-pressure specific heat from pressure and temperature
ADReal
HEM::cp_liquid_from_p_T(ADReal p, ADReal T) const
{
  return _fp_liquid.cp_from_p_T(p, T);
}

// Vapor phase constant-pressure specific heat from pressure and temperature
ADReal
HEM::cp_vapor_from_p_T(ADReal p, ADReal T) const
{
  return _fp_vapor.cp_from_p_T(p, T);
}

// Mixture constant-volume specific heat from pressure and temperature
ADReal
HEM::cv_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const
{
  return (1 - alpha) * (_fp_liquid.cv_from_p_T(p, T)) + alpha * (_fp_vapor.cv_from_p_T(p, T));
}

// Liquid phase constant-volume specific heat from pressure and temperature
ADReal
HEM::cv_liquid_from_p_T(ADReal p, ADReal T) const
{
  return _fp_liquid.cv_from_p_T(p, T);
}

// Vapor phase constant-volume specific heat from pressure and temperature
ADReal
HEM::cv_vapor_from_p_T(ADReal p, ADReal T) const
{
  return _fp_vapor.cv_from_p_T(p, T);
}

// Mixture speed of sound as a function of the pressure and temperature.
ADReal
HEM::c_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const
{
  // defining the necessary liquid properties
  ADReal v_liquid, e_liquid, c_liquid, rho_liquid;

  // defining the necessary vapor properties
  ADReal v_vapor, e_vapor, c_vapor, rho_vapor;

  // defining the Wood's equation for the mixture speed of sound terms
  ADReal C1, C2, C3, C4;

  // calculating v and e from p and T
  _fp_liquid.v_e_from_p_T(p, T, v_liquid, e_liquid);
  _fp_vapor.v_e_from_p_T(p, T, v_vapor, e_vapor);

  // calculating c from v and e
  c_liquid = _fp_liquid.c_from_v_e(v_liquid, e_liquid);
  c_vapor = _fp_vapor.c_from_v_e(v_vapor, e_vapor);

  // calculating fluid and vapor densitites as a function of p and T
  rho_liquid = _fp_liquid.rho_from_p_T(p, T);
  rho_vapor = _fp_vapor.rho_from_p_T(p, T);

  // Calculating the terms for the Wood's equation
  C1 = std::sqrt(alpha * rho_vapor + (1 - alpha) * rho_liquid);
  C2 = rho_liquid * rho_vapor;
  C3 = alpha * rho_liquid * std::pow(c_liquid, 2) + (1 - alpha) * rho_vapor * std::pow(c_vapor, 2);
  C4 = c_vapor * c_liquid;

  // returning the mixture speed of sound
  return (1 / C1) * std::sqrt(C2 / C3) * C4;
}

// Liquid phase speed of sound as a function of the pressure and temperature.
ADReal
HEM::c_liquid_from_p_T(ADReal p, ADReal T) const
{
  // defining the necessary liquid properties
  ADReal v_liquid, e_liquid;

  // calculating v and e from p and T
  _fp_liquid.v_e_from_p_T(p, T, v_liquid, e_liquid);

  // returning the liquid phase speed of sound
  return _fp_liquid.c_from_v_e(v_liquid, e_liquid);
}

// Vapor phase speed of sound as a function of the pressure and temperature.
ADReal
HEM::c_vapor_from_p_T(ADReal p, ADReal T) const
{
  // defining the necessary vapor properties
  ADReal v_vapor, e_vapor;

  // calculating v and e from p and T
  _fp_vapor.v_e_from_p_T(p, T, v_vapor, e_vapor);

  // returning the vapor phase speed of sound
  return _fp_vapor.c_from_v_e(v_vapor, e_vapor);
}

// Liquid phase pressure as a function of the specific volume and internal energy.
ADReal
HEM::p_liquid_from_v_e(ADReal v, ADReal e) const
{
  return _fp_liquid.p_from_v_e(v, e);
}

// Vapor phase pressure as a function of the specific volume and internal energy.
ADReal
HEM::p_vapor_from_v_e(ADReal v, ADReal e) const
{
  return _fp_vapor.p_from_v_e(v, e);
}

// Liquid phase temperature as a function of the specific volume and internal energy.
ADReal
HEM::T_liquid_from_v_e(ADReal v, ADReal e) const
{
  return _fp_liquid.T_from_v_e(v, e);
}

// Vapor phase temperature as a function of the specific volume and internal energy.
ADReal
HEM::T_vapor_from_v_e(ADReal v, ADReal e) const
{
  return _fp_vapor.T_from_v_e(v, e);
}
