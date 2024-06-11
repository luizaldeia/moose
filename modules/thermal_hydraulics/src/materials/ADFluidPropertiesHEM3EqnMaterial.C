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
  HEM::HEMState state = _fp.fluid_state(_rhoA[_qp], _rhouA[_qp], _rhoEA[_qp], _area[_qp]);

  _rho[_qp] = state.rho;
  _v[_qp] = state.v;
  _vel[_qp] = state.vel;
  _e[_qp] = state.e;
  _p[_qp] = state.p;
  _T[_qp] = state.T;
  _h[_qp] = state.h;
  _H[_qp] = state.H;
  _c[_qp] = state.c;
  _cp[_qp] = state.cp;
  _cv[_qp] = state.cv;
  _k[_qp] = state.k;
  _alpha[_qp] = state.alpha;
  // std::cout << "qp ---------------> " << _qp << "\n";
  // std::cout << "\n";
  // std::cout << "========================================\n";
  // std::cout << "ADFluidPropertiesHEM3EqnMaterial \n";
  // std::cout << "========================================\n";
  // std::cout << "\n";
}
