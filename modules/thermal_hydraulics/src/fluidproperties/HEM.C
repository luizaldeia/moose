//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HEM.h"

registerMooseObject("ThermalHydraulicsApp", HEM);

InputParameters
HEM::validParams()
{
  InputParameters params = TwoPhaseFluidProperties::validParams();
  params.addClassDescription("Class for the Homogeneous Equilibrium Model fluid properties as a "
                             "function of pressure and temperature for water");
  params.addRequiredParam<UserObjectName>("fp_2phase", "Two-phase water property user object name");
  params.addParam<Real>("temp_guess", 300, "Initial temperature guess");
  // This is necessary because initialize() must be called before any interface
  // can be used (which can occur as early as initialization of variables).
  params.set<ExecFlagEnum>("execute_on") = EXEC_INITIAL;

  return params;
}

HEM::HEM(const InputParameters & parameters)
  : TwoPhaseFluidProperties(parameters),

    _fp_2phase(getUserObject<TwoPhaseFluidProperties>("fp_2phase")),

    _fp_liquid(getUserObjectByName<SinglePhaseFluidProperties>(_fp_2phase.getLiquidName())),

    _fp_vapor(getUserObjectByName<SinglePhaseFluidProperties>(_fp_2phase.getVaporName())),

    _temp_guess(getParam<Real>("temp_guess"))
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
  ADReal el = (1 - alpha) * (_fp_liquid.e_from_p_T(p, T)) * (_fp_liquid.rho_from_p_T(p, T));
  ADReal eg = alpha * (_fp_vapor.e_from_p_T(p, T)) * (_fp_vapor.rho_from_p_T(p, T));
  return (el + eg) / rho_mixture_from_p_T(p, T, alpha);
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

// Vapor phase internal energy from pressure and density
ADReal
HEM::e_vapor_from_p_rho(ADReal p, ADReal rho) const
{
  return _fp_vapor.e_from_p_rho(p, rho);
}

// Liquid phase internal energy from pressure and density
ADReal
HEM::e_liquid_from_p_rho(ADReal p, ADReal rho) const
{
  return _fp_liquid.e_from_p_rho(p, rho);
}

// Liquid phase internal energy from pressure and density
ADReal
HEM::e_mixture_from_p_rho(ADReal p, ADReal rho, ADReal alpha) const
{
  if (alpha <= 0.0)
  {
    return e_liquid_from_p_rho(p, rho);
  }
  else if (alpha >= 1.0)
  {
    return e_vapor_from_p_rho(p, rho);
  }
  else
  {
    ADReal T = _fp_2phase.T_sat(p);
    return e_mixture_from_p_T(p, T, alpha);
  }
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

// Liquid phase speed of sound as a function of the specific volume and internal energy.
ADReal
HEM::c_liquid_from_v_e(ADReal v, ADReal e) const
{
  // returning the liquid phase speed of sound
  return _fp_liquid.c_from_v_e(v, e);
}

// Vapor phase speed of sound as a function of the specific volume and internal energy.
ADReal
HEM::c_vapor_from_v_e(ADReal v, ADReal e) const
{
  // returning the vapor phase speed of sound
  return _fp_vapor.c_from_v_e(v, e);
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

// Mixture dynamic vicosity from specific volume and internal energy
ADReal
HEM::mu_mixture_from_p_T(ADReal p, ADReal T, ADReal alpha) const
{
  return (1 - alpha) * (_fp_liquid.mu_from_p_T(p, T)) + alpha * (_fp_vapor.mu_from_p_T(p, T));
}

// Liquid phase dynamic vicosity from specific volume and internal energy
ADReal
HEM::mu_liquid_from_v_e(ADReal v, ADReal e) const
{
  return _fp_liquid.mu_from_v_e(v, e);
}

// Vapor phase dynamic vicosity from specific volume and internal energy
ADReal
HEM::mu_vapor_from_v_e(ADReal v, ADReal e) const
{
  return _fp_vapor.mu_from_v_e(v, e);
}

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Calculating the relevant fluid properties as a function of the conservative variables:
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// density time internal energy from rho*E*A, rho*u*A, rho*A, and A: e = rhoe = rhoEA/A -
// 0.5*rho*vel^2
ADReal
HEM::rhoe_from_rhoEA_rhouA_rhoA_A(ADReal rhoEA, ADReal rhouA, ADReal rhoA, ADReal A) const
{
  return rhoEA / A - 0.5 * THM::rho_from_arhoA_alpha_A(rhoA, 1.0, A) *
                         THM::vel_from_arhoA_arhouA(rhouA, rhoA) *
                         THM::vel_from_arhoA_arhouA(rhouA, rhoA);
}

HEM::HEMState
HEM::fluid_state(ADReal rhoA, ADReal rhouA, ADReal rhoEA, ADReal A) const
{
  HEM::HEMState state;

  std::cout << "\n";
  std::cout << "========================================\n";
  std::cout << "State \n";
  std::cout << "rhoA ----------> " << rhoA.value() << "\n";
  std::cout << "rhouA ----------> " << rhouA.value() << "\n";
  std::cout << "rhoEA ----------> " << rhoEA.value() << "\n";
  std::cout << "A --------------> " << A.value() << "\n";
  std::cout << "\n";

  // Declaring some relevant fluid properties
  ADReal rho, v, vel, e, p, T, h, H, c, cp, cv, k, alpha;

  // density
  rho = THM::rho_from_arhoA_alpha_A(rhoA, 1.0, A);
  std::cout << "rho -----------> " << rho.value() << "\n";
  // specific volume
  v = 1.0 / THM::rho_from_arhoA_alpha_A(rhoA, 1.0, A);
  std::cout << "v -------------> " << v.value() << "\n";
  // velocity
  vel = THM::vel_from_arhoA_arhouA(rhoA, rhouA);
  std::cout << "vel -----------> " << vel.value() << "\n";
  // internal energy
  e = THM::e_from_arhoA_arhouA_arhoEA(rhoA, rhouA, rhoEA);
  std::cout << "e -------------> " << e.value() << "\n";

  // Initiating the tolerances and error function
  double tol_e = 1E-4;

  // Initiating the e error:
  ADReal error = 0.0;

  // Newton method to minimize the error function
  ADReal de_dT_old, de_dT, error_de_dT;
  ADReal T_plus, T_minus, p_plus, p_minus;
  ADReal rho_l_plus, rho_g_plus, e_l_plus, e_g_plus, f_plus;
  ADReal rho_l_minus, rho_g_minus, e_l_minus, e_g_minus, f_minus;
  Real epsilon, tol_der;

  tol_der = 1e-4;

  // Declaring some relevant fluid properties to estimate the fluid state (1 or 2-phase)
  ADReal Psat, rho_sat_l, rho_sat_g, e_sat_l, e_sat_g, T_guess, rhoe;

  // Declaring some mixture fluid properties used to estimate the fluid state
  ADReal v_m, e_m, p_m, T_m;

  // Defining the number of iterations counter
  int it = 0;

  // Declaring the intital guess for the fluid temperature
  T_guess = _temp_guess;
  rhoe = rho * e;
  std::cout << "T_guess -------> " << T_guess.value() << "\n";
  do
  {

    if (T_guess < 273.25)
    {
      T_guess = 273.25;
    }
    else if (T_guess > 1274.8)
    {
      T_guess = 1274.8;
    }
    std::cout << "T_guess -------> " << T_guess.value() << "\n";
    std::cout << "p_m -------> " << p_m.value() << "\n";

    // Calculating the saturation pressure at T_guess
    Psat = p_sat(T_guess);

    // Calculating the saturation density for each phase
    rho_sat_l = rho_liquid_from_p_T(Psat, T_guess);
    rho_sat_g = rho_vapor_from_p_T(Psat, T_guess);
    e_sat_l = e_liquid_from_p_T(Psat, T_guess);
    e_sat_g = e_vapor_from_p_T(Psat, T_guess);

    // Calculating the void fraction
    alpha = (rho - rho_sat_l) / (rho_sat_g - rho_sat_l);
    // std::cout << "alpha: " << alpha.value() << "\n";
    if (alpha <= 0)
    {
      alpha = 0.0;
    }
    else if (alpha >= 1.0)
    {
      alpha = 1.0;
    }
    if (T_guess > 646.9)
    {
      alpha = 1.0;
    }

    // The fluid properties will be calculated in different ways according to the fluid phase
    // (Liquid, Vapor, Mixture), for simplicity, the single phase fluid properties bellow were
    // called using the subscript _m.

    if (alpha == 0)
    {
      // if alpha <= 0, we are in a subcooled single-phase flow regime were only the liquid
      // exists.
      std::cout << "========================================\n";
      std::cout << "liquid\n";
      // Pressure from _v and _e
      p_m = p_liquid_from_v_e(v, e);
      std::cout << "p_m --------------> " << p_m.value() << "\n";

      if (std::isnan(p_m) or (p_m < 611.210432))
      {
        T_guess += error;
        p_m = p_sat(T_guess);
      }

      // Temperature from _v and _e
      T_m = T_guess;
      std::cout << "T_m --------------> " << T_m.value() << "\n";

      if (T_m > T_sat(p_m))
      {
        T_m = T_sat(p_m) - 0.1;
      }

      // Calculating the single phase liquid specific volume as a function of p_m and T_m.
      v_m = v_liquid_from_p_T(p_m, T_m);

      // Calculating the single phase liquid internal energy as a function of p_m and T_m.
      e_m = e_liquid_from_p_T(p_m, T_m);

      // solving the error function derivative as function of T_guess

      epsilon = 0.1;
      de_dT_old = 1000.0;
      error_de_dT = 0.0;

      do
      {
        T_plus = T_m + epsilon;
        T_minus = T_m - epsilon;

        e_l_plus = e_liquid_from_p_T(p_m, T_plus);
        e_l_minus = e_liquid_from_p_T(p_m, T_minus);

        f_plus = std::abs((e_l_plus - e) / e);
        f_minus = std::abs((e_l_minus - e) / e);

        de_dT = (f_plus - f_minus) / (2 * epsilon);

        error_de_dT = std::abs(de_dT_old - de_dT);
        de_dT_old = de_dT;
        epsilon *= 0.1;

      } while (error_de_dT > tol_der);
      if (de_dT == 0)
      {
        T_guess = T_liquid_from_v_e(v, e);
      }
      else
      {
        T_guess = T_m - 0.5 * error / de_dT;
      }

      if (T_guess > T_sat(p_m))
      {
        T_guess = T_sat(p_m);
      }
      else if (T_guess < 273.25)
      {
        T_guess = T_sat(p_m);
      }

      error = std::abs((e_m - e) / e);
      std::cout << "p_m --------------> " << p_m.value() << "\n";
      std::cout << "T_guess ----------> " << T_guess.value() << "\n";
      std::cout << "error ----------> " << error.value() << "\n";
      std::cout << "de_dT ------------> " << de_dT.value() << "\n";
    }
    else if (alpha == 1)
    {
      // if alpha >= 0, we are in a superheated single-phase flow regime were only the vapor
      // exists.
      std::cout << "========================================\n";
      std::cout << "vapor\n";
      // Pressure from _v and _e
      p_m = p_vapor_from_v_e(v, e);
      std::cout << "p_m ----------> " << p_m.value() << "\n";
      if ((std::isnan(p_m)) or (p_m < 611.210432))
      {
        p_m = p_sat(T_guess);
        std::cout << "p_m nan ----------> " << p_m.value() << "\n";
      }

      // Temperature
      T_m = T_guess;

      if (T_m < T_sat(p_m))
      {
        T_m = T_sat(p_m);
      }

      // Calculating the single phase vapor specific volume as a function of p_m and T_m.
      v_m = v_vapor_from_p_T(p_m, T_m);

      // Calculating the single phase vapor internal energy as a function of p_m and T_m.
      e_m = e_vapor_from_p_T(p_m, T_m);

      // solving the error function derivative as function of T_guess

      epsilon = 0.1;
      de_dT_old = 1000.0;
      error_de_dT = 0.0;

      do
      {
        T_plus = T_m + epsilon;
        T_minus = T_m - epsilon;

        e_l_plus = e_vapor_from_p_T(p_m, T_plus);
        e_l_minus = e_vapor_from_p_T(p_m, T_minus);

        f_plus = std::abs((e_l_plus - e) / e);
        f_minus = std::abs((e_l_minus - e) / e);

        de_dT = (f_plus - f_minus) / (2 * epsilon);

        error_de_dT = std::abs(de_dT_old - de_dT);
        de_dT_old = de_dT;
        epsilon *= 0.1;

      } while (error_de_dT > tol_der);

      if (de_dT == 0)
      {
        T_guess = T_vapor_from_v_e(v, e);
      }
      else
      {
        T_guess = T_m - 0.5 * error / de_dT;
      }

      if (T_guess < T_sat(p_m))
      {
        T_guess = T_sat(p_m);
      }

      error = std::abs((e_m - e) / e);

      std::cout << "T_guess ----------> " << T_guess.value() << "\n";
      std::cout << "error ----------> " << error.value() << "\n";
      std::cout << "de_dT ------------> " << de_dT.value() << "\n";
    }
    else
    {
      // if 0<alpha<1, we are in a two-phase flow regime were we assumed we have an homogeneous
      // mixture between the two individual phases
      std::cout << "========================================\n";
      std::cout << "mixture\n";

      if (T_guess > 646.9)
      {
        T_guess *= 0.5;
      }

      std::cout << "p_mb -------------> " << p_m.value() << "\n";
      std::cout << "T_mb -------------> " << T_m.value() << "\n";

      // The mixture will be at Psat
      p_m = Psat;

      // The mixture will be at T_guess
      T_m = T_guess;
      if (T_guess > 646.9)
      {
        T_guess *= 0.5;
      }

      std::cout << "p_ma -------------> " << p_m.value() << "\n";
      std::cout << "T_ma -------------> " << T_m.value() << "\n";

      // std::cout << "T_guess ----------> " << T_guess.value() << "\n";
      // Calculating the mixture specific volume as a function of Psat and T_guest.
      v_m = v_mixture_from_p_T(p_m, T_m, alpha);

      // Calculating the mixture internal energy as a function of Psat and T_guest.
      e_m = e_mixture_from_p_T(p_m, T_m, alpha);

      // Calculating the new mixture temperature (for the mixture both phases, liquid and vapor,
      // should have the same temperature).

      // solving the error function derivative as function of T_guess

      epsilon = 0.1;
      de_dT_old = 1000.0;
      error_de_dT = 0.0;

      do
      {
        T_plus = T_guess + epsilon;
        T_minus = T_guess - epsilon;
        p_plus = p_sat(T_plus);
        p_minus = p_sat(T_minus);

        rho_l_plus = rho_liquid_from_p_T(p_plus, T_plus);
        rho_g_plus = rho_vapor_from_p_T(p_plus, T_plus);
        e_l_plus = e_liquid_from_p_T(p_plus, T_plus);
        e_g_plus = e_vapor_from_p_T(p_plus, T_plus);

        ADReal C1_plus =
            (rhoe - rho_l_plus * e_l_plus) / (rho_g_plus * e_g_plus - rho_l_plus * e_l_plus);
        ADReal C2_plus = (rho - rho_l_plus) / (rho_g_plus - rho_l_plus);
        f_plus = std::abs(C1_plus - C2_plus);

        rho_l_minus = rho_liquid_from_p_T(p_minus, T_minus);
        rho_g_minus = rho_vapor_from_p_T(p_minus, T_minus);
        e_l_minus = e_liquid_from_p_T(p_minus, T_minus);
        e_g_minus = e_vapor_from_p_T(p_minus, T_minus);

        ADReal C1_minus =
            (rhoe - rho_l_minus * e_l_minus) / (rho_g_minus * e_g_minus - rho_l_minus * e_l_minus);
        ADReal C2_minus = (rho - rho_l_minus) / (rho_g_minus - rho_l_minus);
        f_minus = std::abs(C1_minus - C2_minus);

        de_dT = (f_plus - f_minus) / (2 * epsilon);
        error_de_dT = std::abs(de_dT_old - de_dT);
        de_dT_old = de_dT;
        epsilon *= 0.1;

      } while (error_de_dT > tol_der);

      T_guess = T_m - 0.5 * error / de_dT;
      std::cout << "T_guess ----------> " << T_guess.value() << "\n";
      if (T_guess > 646.9)
      {
        T_guess *= 0.5;
      }
      else if (T_guess < 273.25)
      {
        T_guess = 273.25;
      }

      ADReal C1 = (rhoe - rho_sat_l * e_sat_l) / (rho_sat_g * e_sat_g - rho_sat_l * e_sat_l);
      ADReal C2 = (rho - rho_sat_l) / (rho_sat_g - rho_sat_l);
      error = std::abs(C1 - C2);

      if (error > tol_e)
      {
        p_m = p_sat(T_guess);
      }

      std::cout << "T_guess ----------> " << T_guess.value() << "\n";
      std::cout << "p_m --------------> " << p_m.value() << "\n";
      std::cout << "error ------------> " << error.value() << "\n";
      std::cout << "de_dT ------------> " << de_dT.value() << "\n";
    }
    if (it >= 100)
    {
      mooseError(
          "The maximum number of iterations was reached while trying to determine the fluid state");
      break;
    }
    it++;
    std::cout << "it ---------------> " << it << "\n";
  } while (error > tol_e);

  // Fluid pressure
  p = p_m;
  std::cout << "p_m --------------> " << p_m.value() << "\n";
  // Fluid temperature
  T = T_m;

  // Fluid specific enthalpy
  h = e + p / rho;

  // Fluid specific total (stagnation) enthalpy
  H = h + 0.5 * vel * vel;

  if (alpha <= 0)
  {
    // Fluid speed of sound
    c = c_liquid_from_p_T(p, T);

    // Fluid constant-pressure specific heat
    cp = cp_liquid_from_p_T(p, T);

    // Fluid constant-volume specific heat
    cv = cv_liquid_from_p_T(p, T);

    // Fluid thermal conductivity
    k = k_liquid_from_p_T(p, T);
  }
  else if (alpha >= 1)
  {
    // Fluid speed of sound
    c = c_vapor_from_p_T(p, T);

    // Fluid constant-pressure specific heat
    cp = cp_vapor_from_p_T(p, T);

    // Fluid constant-volume specific heat
    cv = cv_vapor_from_p_T(p, T);

    // Fluid thermal conductivity
    k = k_vapor_from_p_T(p, T);
  }
  else
  {
    // Fluid speed of sound
    c = c_mixture_from_p_T(p, T, alpha);

    // Fluid constant-pressure specific heat
    cp = cp_mixture_from_p_T(p, T, alpha);

    // Fluid constant-volume specific heat
    cv = cv_mixture_from_p_T(p, T, alpha);

    // Fluid thermal conductivity
    k = k_mixture_from_p_T(p, T, alpha);
  }

  state.rho = rho;
  state.v = v;
  state.vel = vel;
  state.e = e;
  state.p = p;
  state.T = T;
  state.h = h;
  state.H = H;
  state.c = c;
  state.cp = cp;
  state.cv = cv;
  state.k = k;
  state.alpha = alpha;

  return state;
}
