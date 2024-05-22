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

    _v_m(declareADProperty<Real>("v_m")),

    _vel(declareADProperty<Real>("vel")),

    _e(declareADProperty<Real>("e")),

    _e_m(declareADProperty<Real>("e_m")),

    _p(declareADProperty<Real>("p")),

    _T(declareADProperty<Real>("T")),

    _h(declareADProperty<Real>("h")),

    _H(declareADProperty<Real>("H")),

    _c(declareADProperty<Real>("c")),

    _cp(declareADProperty<Real>("cp")),

    _cv(declareADProperty<Real>("cv")),

    _k(declareADProperty<Real>("k")),

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
  _e[_qp] = (_rhoEA[_qp] - 0.5 * _rhouA[_qp] * _rhouA[_qp] / _rhoA[_qp]) / _rhoA[_qp];

  // Initiating the tolerances for v and e
  double tol_v = 1E-6;
  double tol_e = 1E-3;

  // Initiating the errors for v and e.
  ADReal error_v, error_e;

  // Declaring the relevant fluid properties
  ADReal Psat, rho_lsat, rho_vsat, alpha;

  // Defining the number of iterations counter
  int it = 0;
  do
  {
    // Declaring the intital guess for the fluid temperature
    _T[_qp] = 400;

    // Calculating the saturation pressure at T0
    Psat = _fp.p_sat(_T[_qp]);

    // Calculating the saturation density for each phase
    rho_lsat = _fp.rho_liquid_from_p_T(Psat, _T[_qp]);
    rho_vsat = _fp.rho_vapor_from_p_T(Psat, _T[_qp]);

    // Calculating the void fraction
    alpha = (_rho[_qp] - rho_lsat) / (rho_vsat - rho_lsat);

    // The fluid properties will be calculated in different ways according to the fluid phase
    // (Liquid, Vapor, Mixture) and regime (Single- or Two-Phase flow)

    if (alpha <= 0)
    {
      // if alpha <= 0, we are in a subcooled single-phase flow regime were only the liquid exists.

      // Pressure from v_THM and e_THM
      _p[_qp] = _fp.p_liquid_from_v_e(_v[_qp], _e[_qp]);

      // Temperature from v_THM and e_THM
      _T[_qp] = _fp.T_liquid_from_v_e(_v[_qp], _e[_qp]);

      // Calculating the single phase liquid specific volume as a function of p_THM and T_THM. Here
      // it is being caleed as v_m (mixture specific volume) for simplicity.
      _v_m[_qp] = _fp.v_liquid_from_p_T(_p[_qp], _T[_qp]);

      // Calculating the single phase liquid internal energy as a function of p_THM and T_THM. Here
      // it is being caleed as e_m (mixture internal energy) for simplicity.
      _e_m[_qp] = _fp.e_liquid_from_p_T(_p[_qp], _T[_qp]);

      _T[_qp] = _fp.T_liquid_from_v_e(_v_m[_qp], _e_m[_qp]);
    }
    else if (alpha >= 0)
    {
      // if alpha >= 0, we are in a superheated single-phase flow regime were only the vapor exists.

      // Pressure from v_THM and e_THM
      _p[_qp] = _fp.p_vapor_from_v_e(_v[_qp], _e[_qp]);

      // Temperature from v_THM and e_THM
      _T[_qp] = _fp.T_vapor_from_v_e(_v[_qp], _e[_qp]);

      // Calculating the single phase vapor specific volume as a function of p_THM and T_THM. Here
      // it is being caleed as v_m (mixture specific volume) for simplicity.
      _v_m[_qp] = _fp.v_vapor_from_p_T(_p[_qp], _T[_qp]);

      // Calculating the single phase vapor internal energy as a function of p_THM and T_THM. Here
      // it is being caleed as e_m (mixture internal energy) for simplicity.
      _e_m[_qp] = _fp.e_vapor_from_p_T(_p[_qp], _T[_qp]);

      _T[_qp] = _fp.T_vapor_from_v_e(_v_m[_qp], _e_m[_qp]);
    }
    else
    {
      // if 0<alpha<1, we are in a two-phase flow regime were we assumed we have an homogeneous
      // mixture between the two individual phases

      // Calculating the mixture specific volume as a function of Psat and T.
      _v_m[_qp] = _fp.v_mixture_from_p_T(Psat, _T[_qp], alpha);

      // Calculating the mixture internal energy as a function of Psat and T.
      _e_m[_qp] = _fp.e_mixture_from_p_T(Psat, _T[_qp], alpha);

      // Calculating the new mixture temperature.
      _T[_qp] = _fp.T_liquid_from_v_e(_e_m[_qp], _v_m[_qp]);
    }
    if (it >= 1000)
    {
      mooseWarning("The maximum number of iterations was reached while trying to determine the "
                   "fluid state");
    }

    error_v = std::abs(_v[_qp] - _v_m[_qp]);
    error_e = std::abs(_e[_qp] - _e_m[_qp]);
    it++;

  } while (error_v > tol_v && error_e > tol_e);

  _h[_qp] = _e[_qp] + _p[_qp] / _rho[_qp];

  _H[_qp] = _h[_qp] + 0.5 * _vel[_qp] * _vel[_qp];

  _c[_qp] = _fp.c_mixture_from_p_T(_p[_qp], _T[_qp], alpha);

  _cp[_qp] = _fp.cp_mixture_from_p_T(_p[_qp], _T[_qp], alpha);

  _cv[_qp] = _fp.cv_mixture_from_p_T(_p[_qp], _T[_qp], alpha);

  _k[_qp] = _fp.k_mixture_from_p_T(_p[_qp], _T[_qp], alpha);
}
