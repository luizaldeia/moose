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

  // std::cout << "qp --------> " << _qp << "\n";

  // density from rho*A
  _rho[_qp] = _rhoA[_qp] / _area[_qp];

  // specific volume from rho*A
  _v[_qp] = 1.0 / _rho[_qp];

  // flow velocity from rho*u*A and rho*A
  _vel[_qp] = _rhouA[_qp] / _rhoA[_qp];

  // internal energy from rho*E*A, rho*u*A, and rho*A

  ADReal rhoe = _rhoEA[_qp] / _area[_qp] - 0.5 * _rho[_qp] * _vel[_qp] * _vel[_qp];

  _e[_qp] = rhoe / _rho[_qp];

  // Initiating the tolerances for v and e
  // double tol_v = 1E-4;
  double tol_e = 1E-4;

  // Initiating the error from the previous iteration for e.
  // ADReal error_e_old = std::numeric_limits<double>::max();

  // Initiating the e error:
  ADReal error_e = 0.0;

  // Declaring some relevant fluid properties to estimate the fluid state (1 or 2-phase)
  ADReal Psat, rho_sat_l, rho_sat_g, e_sat_l, e_sat_g, alpha, T_guess;

  // Declaring some mixture fluid properties used to estimate the fluid state
  ADReal v_m, e_m, p_m, T_m;

  // Defining the number of iterations counter
  int it = 0;

  // Declaring the intital guess for the fluid temperature
  T_guess = 300;

  // Declaring the intital step size.
  // double step_size = 1.0;

  // Declaring the annealing factor ro reduce the step size
  // double annealing_factor = 0.99;

  do
  {
    // std::cout << "T_guess ---> " << T_guess.value() << "\n";
    std::cout << "it --------> " << it << "\n";

    // Calculating the saturation pressure at T_guess
    Psat = _fp.p_sat(T_guess);
    // std::cout << "Psat-------> " << Psat.value() << "\n";
    // Calculating the saturation density for each phase
    rho_sat_l = _fp.rho_liquid_from_p_T(Psat, T_guess);
    rho_sat_g = _fp.rho_vapor_from_p_T(Psat, T_guess);
    e_sat_l = _fp.e_liquid_from_p_T(Psat, T_guess);
    e_sat_g = _fp.e_vapor_from_p_T(Psat, T_guess);

    // Calculating the void fraction
    alpha = (_rho[_qp] - rho_sat_l) / (rho_sat_g - rho_sat_l);
    // std::cout << "alpha -----> " << alpha.value() << "\n";

    if (alpha < 1E-2)
    {
      alpha = 0.0;
    }
    else if (alpha > 0.99)
    {
      alpha = 1.0;
    }

    // The fluid properties will be calculated in different ways according to the fluid phase
    // (Liquid, Vapor, Mixture), for simplicity, the single phase fluid properties bellow were
    // called using the subscript _m.

    if (alpha == 0)
    {
      // if alpha <= 0, we are in a subcooled single-phase flow regime were only the liquid exists.

      std::cout << "======================================================== \n";
      std::cout << "Liquid \n";

      // Pressure from _v and _e
      p_m = _fp.p_liquid_from_v_e(_v[_qp], _e[_qp]);

      // Temperature from _v and _e
      T_m = T_guess;

      // Calculating the single phase liquid specific volume as a function of p_m and T_m.
      v_m = _fp.v_liquid_from_p_T(p_m, T_m);

      // Calculating the single phase liquid internal energy as a function of p_m and T_m.
      e_m = _fp.e_liquid_from_p_T(p_m, T_m);

      // Updating T_guess and p_guess (the sum with T_thm is just to converge faster)
      T_guess = 0.5 * (_fp.T_liquid_from_v_e(v_m, e_m) + _fp.T_liquid_from_v_e(_v[_qp], _e[_qp]));

      std::cout << "p_m -------> " << p_m.value() << "\n";
      std::cout << "T_m -------> " << T_m.value() << "\n";
      std::cout << "v_m -------> " << v_m.value() << "\n";
      std::cout << "e_m -------> " << e_m.value() << "\n";
      std::cout << "alpha -----> " << alpha.value() << "\n";
      std::cout << "T_guess ---> " << T_guess.value() << "\n";
      std::cout << "======================================================== \n";
    }
    else if (alpha == 1)
    {
      // if alpha >= 0, we are in a superheated single-phase flow regime were only the vapor exists.

      std::cout << "======================================================== \n";
      std::cout << "Vapor \n";

      // Pressure from _v and _e
      p_m = _fp.p_vapor_from_v_e(_v[_qp], _e[_qp]);

      // Temperature
      T_m = T_guess;

      // Calculating the single phase vapor specific volume as a function of p_m and T_m.
      v_m = _fp.v_vapor_from_p_T(p_m, T_m);

      // Calculating the single phase vapor internal energy as a function of p_m and T_m.
      e_m = _fp.e_vapor_from_p_T(p_m, T_m);

      // Updating T_guess and p_guess
      T_guess = 0.5 * (_fp.T_vapor_from_v_e(v_m, e_m) + _fp.T_vapor_from_v_e(_v[_qp], _e[_qp]));

      std::cout << "p_m -------> " << p_m.value() << "\n";
      std::cout << "T_m -------> " << T_m.value() << "\n";
      std::cout << "v_m -------> " << v_m.value() << "\n";
      std::cout << "e_m -------> " << e_m.value() << "\n";
      std::cout << "alpha -----> " << alpha.value() << "\n";
      std::cout << "T_guess ---> " << T_guess.value() << "\n";
      std::cout << "======================================================== \n";
    }
    else
    {
      // if 0<alpha<1, we are in a two-phase flow regime were we assumed we have an homogeneous
      // mixture between the two individual phases

      std::cout << "======================================================== \n";
      std::cout << "Mixture \n";

      // The mixture will be at Psat
      p_m = Psat;

      // The mixture will be at T_guess
      T_m = T_guess;

      // Calculating the mixture specific volume as a function of Psat and T_guest.
      v_m = _fp.v_mixture_from_p_T(p_m, T_m, alpha);

      // Calculating the mixture internal energy as a function of Psat and T_guest.
      e_m = _fp.e_mixture_from_p_T(p_m, T_m, alpha);

      // Calculating the new mixture temperature (for the mixture both phases, liquid and vapor,
      // should have the same temperature).

      // Newton method to minimize the error function
      ADReal de_dT_old, de_dT, error_de_dT;
      // ADReal d2e_dT2_old, d2e_dT2, error_d2e_dT2, g_T, dg_dT;
      ADReal T_sat_plus, T_sat_minus, p_sat_plus, p_sat_minus;
      ADReal rho_sat_l_plus, rho_sat_g_plus, e_sat_l_plus, e_sat_g_plus;
      ADReal rho_sat_l_minus, rho_sat_g_minus, e_sat_l_minus, e_sat_g_minus;
      Real epsilon = 1;
      Real tol_der = 1E-4;
      de_dT_old = 1000;
      // d2e_dT2_old = 1000;
      error_de_dT = 0;
      // error_d2e_dT2 = 0;

      // solving the error function derivative as function of T_guess
      std::cout << "\n";
      std::cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- \n";
      do
      {
        T_sat_plus = T_guess + epsilon;
        T_sat_minus = T_guess - epsilon;
        p_sat_plus = _fp.p_sat(T_sat_plus);
        p_sat_minus = _fp.p_sat(T_sat_minus);

        rho_sat_l_plus = _fp.rho_liquid_from_p_T(p_sat_plus, T_sat_plus);
        rho_sat_g_plus = _fp.rho_vapor_from_p_T(p_sat_plus, T_sat_plus);
        e_sat_l_plus = _fp.e_liquid_from_p_T(p_sat_plus, T_sat_plus);
        e_sat_g_plus = _fp.e_vapor_from_p_T(p_sat_plus, T_sat_plus);

        ADReal C1_plus = (rhoe - rho_sat_l_plus * e_sat_l_plus) /
                         (rho_sat_g_plus * e_sat_g_plus - rho_sat_l_plus * e_sat_l_plus);
        ADReal C2_plus = (_rho[_qp] - rho_sat_l_plus) / (rho_sat_g_plus - rho_sat_l_plus);
        ADReal f_plus = std::abs(C1_plus - C2_plus);

        rho_sat_l_minus = _fp.rho_liquid_from_p_T(p_sat_minus, T_sat_minus);
        rho_sat_g_minus = _fp.rho_vapor_from_p_T(p_sat_minus, T_sat_minus);
        e_sat_l_minus = _fp.e_liquid_from_p_T(p_sat_minus, T_sat_minus);
        e_sat_g_minus = _fp.e_vapor_from_p_T(p_sat_minus, T_sat_minus);

        ADReal C1_minus = (rhoe - rho_sat_l_minus * e_sat_l_minus) /
                          (rho_sat_g_minus * e_sat_g_minus - rho_sat_l_minus * e_sat_l_minus);
        ADReal C2_minus = (_rho[_qp] - rho_sat_l_minus) / (rho_sat_g_minus - rho_sat_l_minus);
        ADReal f_minus = std::abs(C1_minus - C2_minus);

        de_dT = (f_plus - f_minus) / (2 * epsilon);
        // d2e_dT2 = (f_plus + f_minus - 2 * error_e) / (4 * std::pow(epsilon, 2));
        error_de_dT = std::abs(de_dT_old - de_dT);
        // error_d2e_dT2 = std::abs(d2e_dT2_old - d2e_dT2);
        de_dT_old = de_dT;
        // d2e_dT2_old = d2e_dT2;
        epsilon *= 0.1;
        std::cout << "de_dT -----> " << de_dT.value() << "\n";
        std::cout << "e_de_dT ---> " << error_de_dT.value() << "\n";
        // std::cout << "d2e_dT2 -----> " << d2e_dT2.value() << "\n";
        // std::cout << "e_d2e_dT2 ---> " << error_d2e_dT2.value() << "\n";
        // } while ((error_de_dT > tol_der) && (error_d2e_dT2 > tol_der));
      } while (error_de_dT > tol_der);
      // g_T = error_e / de_dT;
      // dg_dT = 1 - (d2e_dT2 / de_dT) * g_T;
      // std::cout << "g_T -------> " << g_T.value() << "\n";
      // std::cout << "dg_dT -----> " << dg_dT.value() << "\n";
      std::cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- \n";
      std::cout << "\n";

      std::cout << "error_e ---> " << error_e.value() << "\n";
      T_guess = T_m - 0.5 * error_e / de_dT;
      // T_guess = T_m - g_T / dg_dT;
      std::cout << "T_guess ---> " << T_guess.value() << "\n";

      if (T_guess > 647)
      {
        T_guess *= 0.5;
      }

      // step_size *= annealing_factor; // Gradually reduce step size

      std::cout << "p_m -------> " << p_m.value() << "\n";
      std::cout << "T_m -------> " << T_m.value() << "\n";
      std::cout << "v_m -------> " << v_m.value() << "\n";
      std::cout << "e_m -------> " << e_m.value() << "\n";
      std::cout << "alpha -----> " << alpha.value() << "\n";
      std::cout << "T_guess ---> " << T_guess.value() << "\n";
      std::cout << "======================================================== \n";
    }
    if (it >= 100)
    {
      mooseError("The maximum number of iterations was reached while trying to determine the "
                 "fluid state");
      break;
    }

    ADReal C1 = (rhoe - rho_sat_l * e_sat_l) / (rho_sat_g * e_sat_g - rho_sat_l * e_sat_l);
    ADReal C2 = (_rho[_qp] - rho_sat_l) / (rho_sat_g - rho_sat_l);
    // error_e_old = error_e;
    // error_e = std::abs((_e[_qp].value() - e_m.value()) / _e[_qp].value());
    error_e = std::abs(C1 - C2);

    it++;
    std::cout << "error_e ---> " << error_e.value() << "\n ";
  } while (error_e > tol_e);

  // Fluid pressure
  _p[_qp] = p_m;
  std::cout << "p ---------> " << _p[_qp].value() << "\n";

  // Fluid temperature
  _T[_qp] = T_m;
  std::cout << "T ---------> " << _T[_qp].value() << "\n";

  // Fluid specific enthalpy
  _h[_qp] = _e[_qp] + _p[_qp] / _rho[_qp];

  // Fluid specific total (stagnation) enthalpy
  _H[_qp] = _h[_qp] + 0.5 * _vel[_qp] * _vel[_qp];

  // Fluid speed of sound
  _c[_qp] = _fp.c_mixture_from_p_T(_p[_qp], _T[_qp], alpha);

  // Fluid constant-pressure specific heat
  _cp[_qp] = _fp.cp_mixture_from_p_T(_p[_qp], _T[_qp], alpha);

  // Fluid constant-volume specific heat
  _cv[_qp] = _fp.cv_mixture_from_p_T(_p[_qp], _T[_qp], alpha);

  // Fluid thermal conductivity
  _k[_qp] = _fp.k_mixture_from_p_T(_p[_qp], _T[_qp], alpha);

  // Fluid void fraction
  _alpha[_qp] = alpha;
}
