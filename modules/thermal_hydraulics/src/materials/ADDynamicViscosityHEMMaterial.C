//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADDynamicViscosityHEMMaterial.h"

registerMooseObject("ThermalHydraulicsApp", ADDynamicViscosityHEMMaterial);

InputParameters
ADDynamicViscosityHEMMaterial::validParams()
{
  InputParameters params = Material::validParams();

  params.addRequiredParam<MaterialPropertyName>("mu", "Dynamic viscosity property");
  params.addRequiredParam<MaterialPropertyName>("v", "Specific volume property");
  params.addRequiredParam<MaterialPropertyName>("e", "Specific internal energy property");
  params.addRequiredParam<MaterialPropertyName>("p", "Pressure");
  params.addRequiredParam<MaterialPropertyName>("T", "Temperature");
  params.addRequiredParam<MaterialPropertyName>("alpha", "Fluid void fraction");

  params.addRequiredParam<UserObjectName>("fp_hem", "HEM fluid properties");

  params.addClassDescription("Computes dynamic viscosity as a material property");

  return params;
}

ADDynamicViscosityHEMMaterial::ADDynamicViscosityHEMMaterial(const InputParameters & parameters)
  : Material(parameters),

    _mu_name(getParam<MaterialPropertyName>("mu")),
    _mu(declareADProperty<Real>(_mu_name)),

    _v(getADMaterialProperty<Real>("v")),

    _e(getADMaterialProperty<Real>("e")),

    _p(getADMaterialProperty<Real>("p")),

    _T(getADMaterialProperty<Real>("T")),

    _alpha(getADMaterialProperty<Real>("alpha")),

    _fp_hem(getUserObject<HEM>("fp_hem"))
{
}

void
ADDynamicViscosityHEMMaterial::computeQpProperties()
{
  // std::cout << "\n";
  // std::cout << "========================================\n";
  // std::cout << "ADDynamicViscosityHEMMaterial \n";
  // std::cout << "========================================\n";
  // std::cout << "\n";
  if (_alpha[_qp] <= 0.0)
  {
    _mu[_qp] = _fp_hem.mu_liquid_from_v_e(_v[_qp], _e[_qp]);
  }
  else if (_alpha[_qp] >= 1.0)
  {
    _mu[_qp] = _fp_hem.mu_vapor_from_v_e(_v[_qp], _e[_qp]);
  }
  else
  {
    _mu[_qp] = _fp_hem.mu_mixture_from_p_T(_p[_qp], _T[_qp], _alpha[_qp]);
  }
}
