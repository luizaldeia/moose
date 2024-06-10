//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADFluidPropertiesHEM3EqnMaterial.h"

registerMooseObject("ThermalHydraulicsApp", ADFluidPropertiesHEM3EqnMaterial);

InputParameters
ADFluidPropertiesHEM3EqnMaterial::validParams()
{
  InputParameters params = Material::validParams();

  params.addRequiredCoupledVar("A", "Cross-sectional area");
  params.addRequiredCoupledVar("rhoA", "Conserved density");
  params.addRequiredCoupledVar("rhouA", "Conserved momentum");
  params.addRequiredCoupledVar("rhoEA", "Conserved total energy");

  params.addRequiredParam<UserObjectName>("fp", "The name of the user object for fluid properties");

  params.addClassDescription(
      "Defines the HEM material properties from fluid properties to serve in the 3-equation model");

  return params;
}

ADFluidPropertiesHEM3EqnMaterial::ADFluidPropertiesHEM3EqnMaterial(
    const InputParameters & parameters)
  : Material(parameters),

    _area(adCoupledValue("A")),

    _rhoA(adCoupledValue("rhoA")),

    _rhouA(adCoupledValue("rhouA")),

    _rhoEA(adCoupledValue("rhoEA")),

    _rho(declareADProperty<Real>("rho")),

    _v(declareADProperty<Real>("v")),

    _vel(declareADProperty<Real>("vel")),

    _e(declareADProperty<Real>("e")),

    _p(declareADProperty<Real>("p")),

    _T(declareADProperty<Real>("T")),

    _h(declareADProperty<Real>("h")),

    _H(declareADProperty<Real>("H")),

    _c(declareADProperty<Real>("c")),

    _cp(declareADProperty<Real>("cp")),

    _cv(declareADProperty<Real>("cv")),

    _k(declareADProperty<Real>("k")),

    _alpha(declareADProperty<Real>("alpha")),

    _fp(getUserObject<HEM>("fp"))
{
}

// Calculate the void fraction and the mixture fluid properties
void
ADFluidPropertiesHEM3EqnMaterial::computeQpProperties()
{
  // density from rho*A
  _rho[_qp] = _rhoA[_qp] / _area[_qp];

  // specific volume from rho*A
  _v[_qp] = 1.0 / _rho[_qp];

  // flow velocity from rho*u*A and rho*A
  _vel[_qp] = _rhouA[_qp] / _rhoA[_qp];

  // internal energy from rho*E*A, rho*u*A, and rho*A

  ADReal rhoe = _rhoEA[_qp] / _area[_qp] - 0.5 * _rho[_qp] * _vel[_qp] * _vel[_qp];

  _e[_qp] = rhoe / _rho[_qp];

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
  ADReal Psat, rho_sat_l, rho_sat_g, e_sat_l, e_sat_g, alpha, T_guess;

  // Declaring some mixture fluid properties used to estimate the fluid state
  ADReal v_m, e_m, p_m, T_m;

  // Defining the number of iterations counter
  int it = 0;

  // if (_qp == 0)
  // {
  // Declaring the intital guess for the fluid temperature
  T_guess = 300.0;

  do
  {

    // Calculating the saturation pressure at T_guess
    Psat = _fp.p_sat(T_guess);

    // Calculating the saturation density for each phase
    rho_sat_l = _fp.rho_liquid_from_p_T(Psat, T_guess);
    rho_sat_g = _fp.rho_vapor_from_p_T(Psat, T_guess);
    e_sat_l = _fp.e_liquid_from_p_T(Psat, T_guess);
    e_sat_g = _fp.e_vapor_from_p_T(Psat, T_guess);

    // Calculating the void fraction
    alpha = (_rho[_qp] - rho_sat_l) / (rho_sat_g - rho_sat_l);
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
      p_m = _fp.p_liquid_from_v_e(_v[_qp], _e[_qp]);
      std::cout << "p_m --------------> " << p_m.value() << "\n";
      if (std::isnan(p_m))
      {
        p_m = _fp.p_sat(T_guess);
      }

      // Temperature from _v and _e
      T_m = T_guess;
      std::cout << "T_m --------------> " << T_m.value() << "\n";

      if (T_m > _fp.T_sat(p_m))
      {
        T_m = _fp.T_sat(p_m) - 0.1;
      }
      else if (T_m < 273.25)
      {
        T_m = 273.25;
      }

      // Calculating the single phase liquid specific volume as a function of p_m and T_m.
      v_m = _fp.v_liquid_from_p_T(p_m, T_m);

      // Calculating the single phase liquid internal energy as a function of p_m and T_m.
      e_m = _fp.e_liquid_from_p_T(p_m, T_m);

      // solving the error function derivative as function of T_guess

      epsilon = 0.1;
      de_dT_old = 1000.0;
      error_de_dT = 0.0;

      do
      {
        T_plus = T_m + epsilon;
        T_minus = T_m - epsilon;

        e_l_plus = _fp.e_liquid_from_p_T(p_m, T_plus);
        e_l_minus = _fp.e_liquid_from_p_T(p_m, T_minus);

        f_plus = std::abs((e_l_plus - _e[_qp]) / _e[_qp]);
        f_minus = std::abs((e_l_minus - _e[_qp]) / _e[_qp]);

        de_dT = (f_plus - f_minus) / (2 * epsilon);

        error_de_dT = std::abs(de_dT_old - de_dT);
        de_dT_old = de_dT;
        epsilon *= 0.1;

      } while (error_de_dT > tol_der);
      if (de_dT == 0)
      {
        T_guess = _fp.T_liquid_from_v_e(_v[_qp], _e[_qp]);
      }
      else
      {
        T_guess = T_m - 0.5 * error / de_dT;
      }

      if (T_guess > _fp.T_sat(p_m))
      {
        T_guess = _fp.T_sat(p_m);
      }

      error = std::abs((e_m - _e[_qp]) / _e[_qp]);

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
      p_m = _fp.p_vapor_from_v_e(_v[_qp], _e[_qp]);

      if (std::isnan(p_m))
      {
        p_m = _fp.p_sat(T_guess);
      }

      // Temperature
      T_m = T_guess;

      if (T_m < _fp.T_sat(p_m))
      {
        T_m = _fp.T_sat(p_m);
      }

      // Calculating the single phase vapor specific volume as a function of p_m and T_m.
      v_m = _fp.v_vapor_from_p_T(p_m, T_m);

      // Calculating the single phase vapor internal energy as a function of p_m and T_m.
      e_m = _fp.e_vapor_from_p_T(p_m, T_m);

      // solving the error function derivative as function of T_guess

      epsilon = 0.1;
      de_dT_old = 1000.0;
      error_de_dT = 0.0;

      do
      {
        T_plus = T_m + epsilon;
        T_minus = T_m - epsilon;

        e_l_plus = _fp.e_vapor_from_p_T(p_m, T_plus);
        e_l_minus = _fp.e_vapor_from_p_T(p_m, T_minus);

        f_plus = std::abs((e_l_plus - _e[_qp]) / _e[_qp]);
        f_minus = std::abs((e_l_minus - _e[_qp]) / _e[_qp]);

        de_dT = (f_plus - f_minus) / (2 * epsilon);

        error_de_dT = std::abs(de_dT_old - de_dT);
        de_dT_old = de_dT;
        epsilon *= 0.1;

      } while (error_de_dT > tol_der);

      if (de_dT == 0)
      {
        T_guess = _fp.T_vapor_from_v_e(_v[_qp], _e[_qp]);
      }
      else
      {
        T_guess = T_m - 0.5 * error / de_dT;
      }

      if (T_guess < _fp.T_sat(p_m))
      {
        T_guess = _fp.T_sat(p_m);
      }
      else if (T_guess > 1274.8)
      {
        T_guess = 1274.8;
      }

      error = std::abs((e_m - _e[_qp]) / _e[_qp]);

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
      // The mixture will be at Psat
      p_m = Psat;

      // The mixture will be at T_guess
      T_m = T_guess;
      if (T_guess > 646.9)
      {
        T_guess = 0.5;
      }

      // std::cout << "T_guess ----------> " << T_guess.value() << "\n";
      // Calculating the mixture specific volume as a function of Psat and T_guest.
      v_m = _fp.v_mixture_from_p_T(p_m, T_m, alpha);

      // Calculating the mixture internal energy as a function of Psat and T_guest.
      e_m = _fp.e_mixture_from_p_T(p_m, T_m, alpha);

      // Calculating the new mixture temperature (for the mixture both phases, liquid and vapor,
      // should have the same temperature).

      // solving the error function derivative as function of T_guess

      epsilon = 1.0;
      de_dT_old = 1000.0;
      error_de_dT = 0.0;

      do
      {
        T_plus = T_guess + epsilon;
        T_minus = T_guess - epsilon;
        p_plus = _fp.p_sat(T_plus);
        p_minus = _fp.p_sat(T_minus);

        rho_l_plus = _fp.rho_liquid_from_p_T(p_plus, T_plus);
        rho_g_plus = _fp.rho_vapor_from_p_T(p_plus, T_plus);
        e_l_plus = _fp.e_liquid_from_p_T(p_plus, T_plus);
        e_g_plus = _fp.e_vapor_from_p_T(p_plus, T_plus);

        ADReal C1_plus =
            (rhoe - rho_l_plus * e_l_plus) / (rho_g_plus * e_g_plus - rho_l_plus * e_l_plus);
        ADReal C2_plus = (_rho[_qp] - rho_l_plus) / (rho_g_plus - rho_l_plus);
        f_plus = std::abs(C1_plus - C2_plus);

        rho_l_minus = _fp.rho_liquid_from_p_T(p_minus, T_minus);
        rho_g_minus = _fp.rho_vapor_from_p_T(p_minus, T_minus);
        e_l_minus = _fp.e_liquid_from_p_T(p_minus, T_minus);
        e_g_minus = _fp.e_vapor_from_p_T(p_minus, T_minus);

        ADReal C1_minus =
            (rhoe - rho_l_minus * e_l_minus) / (rho_g_minus * e_g_minus - rho_l_minus * e_l_minus);
        ADReal C2_minus = (_rho[_qp] - rho_l_minus) / (rho_g_minus - rho_l_minus);
        f_minus = std::abs(C1_minus - C2_minus);

        de_dT = (f_plus - f_minus) / (2 * epsilon);
        error_de_dT = std::abs(de_dT_old - de_dT);
        de_dT_old = de_dT;
        epsilon *= 0.1;

      } while (error_de_dT > tol_der);

      T_guess = T_m - 0.5 * error / de_dT;

      if (T_guess > 646.9)
      {
        T_guess *= 0.5;
      }
      ADReal C1 = (rhoe - rho_sat_l * e_sat_l) / (rho_sat_g * e_sat_g - rho_sat_l * e_sat_l);
      ADReal C2 = (_rho[_qp] - rho_sat_l) / (rho_sat_g - rho_sat_l);
      error = std::abs(C1 - C2);

      std::cout << "T_guess ----------> " << T_guess.value() << "\n";
      std::cout << "error ----------> " << error.value() << "\n";
      std::cout << "de_dT ------------> " << de_dT.value() << "\n";
    }
    if (it >= 100)
    {
      mooseError("The maximum number of iterations was reached while trying to determine the "
                 "fluid state");
      break;
    }
    it++;
    std::cout << "it ---------------> " << it << "\n";
  } while (error > tol_e);

  // Fluid pressure
  _p[_qp] = p_m;

  // Fluid temperature
  _T[_qp] = T_m;

  // Fluid specific enthalpy
  _h[_qp] = _e[_qp] + _p[_qp] / _rho[_qp];

  // Fluid specific total (stagnation) enthalpy
  _H[_qp] = _h[_qp] + 0.5 * _vel[_qp] * _vel[_qp];

  if (alpha <= 0)
  {
    // Fluid speed of sound
    _c[_qp] = _fp.c_liquid_from_p_T(_p[_qp], _T[_qp]);

    // Fluid constant-pressure specific heat
    _cp[_qp] = _fp.cp_liquid_from_p_T(_p[_qp], _T[_qp]);

    // Fluid constant-volume specific heat
    _cv[_qp] = _fp.cv_liquid_from_p_T(_p[_qp], _T[_qp]);

    // Fluid thermal conductivity
    _k[_qp] = _fp.k_liquid_from_p_T(_p[_qp], _T[_qp]);
  }
  else if (alpha >= 1)
  {
    // Fluid speed of sound
    _c[_qp] = _fp.c_vapor_from_p_T(_p[_qp], _T[_qp]);

    // Fluid constant-pressure specific heat
    _cp[_qp] = _fp.cp_vapor_from_p_T(_p[_qp], _T[_qp]);

    // Fluid constant-volume specific heat
    _cv[_qp] = _fp.cv_vapor_from_p_T(_p[_qp], _T[_qp]);

    // Fluid thermal conductivity
    _k[_qp] = _fp.k_vapor_from_p_T(_p[_qp], _T[_qp]);
  }
  else
  {
    // Fluid speed of sound
    _c[_qp] = _fp.c_mixture_from_p_T(_p[_qp], _T[_qp], alpha);

    // Fluid constant-pressure specific heat
    _cp[_qp] = _fp.cp_mixture_from_p_T(_p[_qp], _T[_qp], alpha);

    // Fluid constant-volume specific heat
    _cv[_qp] = _fp.cv_mixture_from_p_T(_p[_qp], _T[_qp], alpha);

    // Fluid thermal conductivity
    _k[_qp] = _fp.k_mixture_from_p_T(_p[_qp], _T[_qp], alpha);
  }

  // Fluid void fraction
  _alpha[_qp] = alpha;
  // }
  // else
  //   // Declaring the intital guess for the fluid temperature
  //   T_guess = _T[_qp - 1];

  // do
  // {

  //   // Calculating the saturation pressure at T_guess
  //   Psat = _fp.p_sat(T_guess);

  //   // Calculating the saturation density for each phase
  //   rho_sat_l = _fp.rho_liquid_from_p_T(Psat, T_guess);
  //   rho_sat_g = _fp.rho_vapor_from_p_T(Psat, T_guess);
  //   e_sat_l = _fp.e_liquid_from_p_T(Psat, T_guess);
  //   e_sat_g = _fp.e_vapor_from_p_T(Psat, T_guess);

  //   // Calculating the void fraction
  //   alpha = (_rho[_qp] - rho_sat_l) / (rho_sat_g - rho_sat_l);
  //   // std::cout << "alpha: " << alpha.value() << "\n";
  //   if (alpha < 5E-3)
  //   {
  //     alpha = 0.0;
  //   }
  //   else if (alpha > 0.995)
  //   {
  //     alpha = 1.0;
  //   }

  //   // The fluid properties will be calculated in different ways according to the fluid phase
  //   // (Liquid, Vapor, Mixture), for simplicity, the single phase fluid properties bellow were
  //   // called using the subscript _m.

  //   if (alpha == 0)
  //   {
  //     // if alpha <= 0, we are in a subcooled single-phase flow regime were only the liquid
  //     // exists.
  //     // std::cout << "========================================\n";
  //     // std::cout << "liquid\n";
  //     // Pressure from _v and _e
  //     p_m = _fp.p_liquid_from_v_e(_v[_qp], _e[_qp]);

  //     // Temperature from _v and _e
  //     T_m = T_guess;

  //     // Calculating the single phase liquid specific volume as a function of p_m and T_m.
  //     v_m = _fp.v_liquid_from_p_T(p_m, T_m);

  //     // Calculating the single phase liquid internal energy as a function of p_m and T_m.
  //     e_m = _fp.e_liquid_from_p_T(p_m, T_m);

  //     // solving the error function derivative as function of T_guess

  //     epsilon = 1.0;
  //     de_dT_old = 1000.0;
  //     error_de_dT = 0.0;

  //     do
  //     {
  //       T_plus = T_m + epsilon;
  //       T_minus = T_m - epsilon;

  //       e_l_plus = _fp.e_liquid_from_p_T(p_m, T_plus);
  //       e_l_minus = _fp.e_liquid_from_p_T(p_m, T_minus);

  //       f_plus = std::abs((e_l_plus - _e[_qp]) / _e[_qp]);
  //       f_minus = std::abs((e_l_minus - _e[_qp]) / _e[_qp]);

  //       de_dT = (f_plus - f_minus) / (2 * epsilon);

  //       error_de_dT = std::abs(de_dT_old - de_dT);
  //       de_dT_old = de_dT;
  //       epsilon *= 0.1;

  //     } while (error_de_dT > tol_der);
  //     if (de_dT == 0)
  //     {
  //       T_guess = _fp.T_liquid_from_v_e(_v[_qp], _e[_qp]);
  //     }
  //     else
  //     {
  //       T_guess = T_m - 0.5 * error / de_dT;
  //     }

  //     if (T_guess > _fp.T_sat(p_m))
  //     {
  //       T_guess = _fp.T_sat(p_m);
  //     }

  //     error = std::abs((e_m - _e[_qp]) / _e[_qp]);

  //     // std::cout << "T_guess ----------> " << T_guess.value() << "\n";
  //     // std::cout << "error ----------> " << error.value() << "\n";
  //     // std::cout << "de_dT ------------> " << de_dT.value() << "\n";
  //   }
  //   else if (alpha == 1)
  //   {
  //     // if alpha >= 0, we are in a superheated single-phase flow regime were only the vapor
  //     // exists.
  //     // std::cout << "========================================\n";
  //     // std::cout << "vapor\n";
  //     // Pressure from _v and _e
  //     p_m = _fp.p_vapor_from_v_e(_v[_qp], _e[_qp]);

  //     // Temperature
  //     T_m = T_guess;

  //     // Calculating the single phase vapor specific volume as a function of p_m and T_m.
  //     v_m = _fp.v_vapor_from_p_T(p_m, T_m);

  //     // Calculating the single phase vapor internal energy as a function of p_m and T_m.
  //     e_m = _fp.e_vapor_from_p_T(p_m, T_m);

  //     // solving the error function derivative as function of T_guess

  //     epsilon = 1.0;
  //     de_dT_old = 1000.0;
  //     error_de_dT = 0.0;

  //     do
  //     {
  //       T_plus = T_m + epsilon;
  //       T_minus = T_m - epsilon;

  //       e_l_plus = _fp.e_vapor_from_p_T(p_m, T_plus);
  //       e_l_minus = _fp.e_vapor_from_p_T(p_m, T_minus);

  //       f_plus = std::abs((e_l_plus - _e[_qp]) / _e[_qp]);
  //       f_minus = std::abs((e_l_minus - _e[_qp]) / _e[_qp]);

  //       de_dT = (f_plus - f_minus) / (2 * epsilon);

  //       error_de_dT = std::abs(de_dT_old - de_dT);
  //       de_dT_old = de_dT;
  //       epsilon *= 0.1;

  //     } while (error_de_dT > tol_der);

  //     if (de_dT == 0)
  //     {
  //       T_guess = _fp.T_vapor_from_v_e(_v[_qp], _e[_qp]);
  //     }
  //     else
  //     {
  //       T_guess = T_m - 0.5 * error / de_dT;
  //     }

  //     if (T_guess < _fp.T_sat(p_m))
  //     {
  //       T_guess = _fp.T_sat(p_m);
  //     }

  //     error = std::abs((e_m - _e[_qp]) / _e[_qp]);

  //     // std::cout << "T_guess ----------> " << T_guess.value() << "\n";
  //     // std::cout << "error ----------> " << error.value() << "\n";
  //     // std::cout << "de_dT ------------> " << de_dT.value() << "\n";
  //   }

  //   else
  //   {
  //     // if 0<alpha<1, we are in a two-phase flow regime were we assumed we have an homogeneous
  //     // mixture between the two individual phases
  //     // std::cout << "========================================\n";
  //     // std::cout << "mixture\n";
  //     // The mixture will be at Psat
  //     p_m = Psat;

  //     // The mixture will be at T_guess
  //     T_m = T_guess;
  //     // std::cout << "T_guess ----------> " << T_guess.value() << "\n";
  //     // Calculating the mixture specific volume as a function of Psat and T_guest.
  //     v_m = _fp.v_mixture_from_p_T(p_m, T_m, alpha);

  //     // Calculating the mixture internal energy as a function of Psat and T_guest.
  //     e_m = _fp.e_mixture_from_p_T(p_m, T_m, alpha);

  //     // Calculating the new mixture temperature (for the mixture both phases, liquid and vapor,
  //     // should have the same temperature).

  //     // solving the error function derivative as function of T_guess

  //     epsilon = 1.0;
  //     de_dT_old = 1000.0;
  //     error_de_dT = 0.0;

  //     do
  //     {
  //       T_plus = T_guess + epsilon;
  //       T_minus = T_guess - epsilon;
  //       p_plus = _fp.p_sat(T_plus);
  //       p_minus = _fp.p_sat(T_minus);

  //       rho_l_plus = _fp.rho_liquid_from_p_T(p_plus, T_plus);
  //       rho_g_plus = _fp.rho_vapor_from_p_T(p_plus, T_plus);
  //       e_l_plus = _fp.e_liquid_from_p_T(p_plus, T_plus);
  //       e_g_plus = _fp.e_vapor_from_p_T(p_plus, T_plus);

  //       ADReal C1_plus =
  //           (rhoe - rho_l_plus * e_l_plus) / (rho_g_plus * e_g_plus - rho_l_plus * e_l_plus);
  //       ADReal C2_plus = (_rho[_qp] - rho_l_plus) / (rho_g_plus - rho_l_plus);
  //       f_plus = std::abs(C1_plus - C2_plus);

  //       rho_l_minus = _fp.rho_liquid_from_p_T(p_minus, T_minus);
  //       rho_g_minus = _fp.rho_vapor_from_p_T(p_minus, T_minus);
  //       e_l_minus = _fp.e_liquid_from_p_T(p_minus, T_minus);
  //       e_g_minus = _fp.e_vapor_from_p_T(p_minus, T_minus);

  //       ADReal C1_minus =
  //           (rhoe - rho_l_minus * e_l_minus) / (rho_g_minus * e_g_minus - rho_l_minus *
  //           e_l_minus);
  //       ADReal C2_minus = (_rho[_qp] - rho_l_minus) / (rho_g_minus - rho_l_minus);
  //       f_minus = std::abs(C1_minus - C2_minus);

  //       de_dT = (f_plus - f_minus) / (2 * epsilon);
  //       error_de_dT = std::abs(de_dT_old - de_dT);
  //       de_dT_old = de_dT;
  //       epsilon *= 0.1;

  //     } while (error_de_dT > tol_der);

  //     T_guess = T_m - 0.5 * error / de_dT;

  //     if (T_guess > 647)
  //     {
  //       T_guess *= 0.5;
  //     }
  //     ADReal C1 = (rhoe - rho_sat_l * e_sat_l) / (rho_sat_g * e_sat_g - rho_sat_l * e_sat_l);
  //     ADReal C2 = (_rho[_qp] - rho_sat_l) / (rho_sat_g - rho_sat_l);
  //     error = std::abs(C1 - C2);

  //     // std::cout << "T_guess ----------> " << T_guess.value() << "\n";
  //     // std::cout << "error ----------> " << error.value() << "\n";
  //     // std::cout << "de_dT ------------> " << de_dT.value() << "\n";
  //   }
  //   if (it >= 100)
  //   {
  //     mooseError("The maximum number of iterations was reached while trying to determine the "
  //                "fluid state");
  //     break;
  //   }
  //   it++;
  //   // std::cout << "it ---------------> " << it << "\n";
  // } while (error > tol_e);

  // // Fluid pressure
  // _p[_qp] = p_m;

  // // Fluid temperature
  // _T[_qp] = T_m;

  // // Fluid specific enthalpy
  // _h[_qp] = _e[_qp] + _p[_qp] / _rho[_qp];

  // // Fluid specific total (stagnation) enthalpy
  // _H[_qp] = _h[_qp] + 0.5 * _vel[_qp] * _vel[_qp];

  // if (_alpha[_qp] <= 0)
  // {
  //   // Fluid speed of sound
  //   _c[_qp] = _fp.c_liquid_from_p_T(_p[_qp], _T[_qp]);

  //   // Fluid constant-pressure specific heat
  //   _cp[_qp] = _fp.cp_liquid_from_p_T(_p[_qp], _T[_qp]);

  //   // Fluid constant-volume specific heat
  //   _cv[_qp] = _fp.cv_liquid_from_p_T(_p[_qp], _T[_qp]);

  //   // Fluid thermal conductivity
  //   _k[_qp] = _fp.k_liquid_from_p_T(_p[_qp], _T[_qp]);
  // }
  // else if (_alpha[_qp] >= 1)
  // {
  //   // Fluid speed of sound
  //   _c[_qp] = _fp.c_vapor_from_p_T(_p[_qp], _T[_qp]);

  //   // Fluid constant-pressure specific heat
  //   _cp[_qp] = _fp.cp_vapor_from_p_T(_p[_qp], _T[_qp]);

  //   // Fluid constant-volume specific heat
  //   _cv[_qp] = _fp.cv_vapor_from_p_T(_p[_qp], _T[_qp]);

  //   // Fluid thermal conductivity
  //   _k[_qp] = _fp.k_vapor_from_p_T(_p[_qp], _T[_qp]);
  // }
  // else
  // {
  //   // Fluid speed of sound
  //   _c[_qp] = _fp.c_mixture_from_p_T(_p[_qp], _T[_qp], alpha);

  //   // Fluid constant-pressure specific heat
  //   _cp[_qp] = _fp.cp_mixture_from_p_T(_p[_qp], _T[_qp], alpha);

  //   // Fluid constant-volume specific heat
  //   _cv[_qp] = _fp.cv_mixture_from_p_T(_p[_qp], _T[_qp], alpha);

  //   // Fluid thermal conductivity
  //   _k[_qp] = _fp.k_mixture_from_p_T(_p[_qp], _T[_qp], alpha);
  // }

  // // Fluid void fraction
  // _alpha[_qp] = alpha;
}
