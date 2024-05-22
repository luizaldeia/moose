// //* This file is part of the MOOSE framework
// //* https://www.mooseframework.org
// //*
// //* All rights reserved, see COPYRIGHT for full restrictions
// //* https://github.com/idaholab/moose/blob/master/COPYRIGHT
// //*
// //* Licensed under LGPL 2.1, please see LICENSE for details
// //* https://www.gnu.org/licenses/lgpl-2.1.html

// #include "ADDynamicViscosityHEMMaterial.h"

// registerMooseObject("ThermalHydraulicsApp", ADDynamicViscosityHEMMaterial);

// InputParameters
// ADDynamicViscosityHEMMaterial::validParams()
// {
//   InputParameters params = Material::validParams();

//   params.addRequiredParam<MaterialPropertyName>("mu", "Dynamic viscosity property");
//   params.addRequiredParam<MaterialPropertyName>("v", "Specific volume property");
//   params.addRequiredParam<MaterialPropertyName>("e", "Specific internal energy property");

//   params.addRequiredParam<UserObjectName>("fp_2phase", "Single-phase fluid properties");

//   params.addClassDescription("Computes dynamic viscosity as a material property");

//   return params;
// }

// ADDynamicViscosityHEMMaterial::ADDynamicViscosityHEMMaterial(const InputParameters & parameters)
//   : Material(parameters),

//     _mu_name(getParam<MaterialPropertyName>("mu")),
//     _mu(declareADProperty<Real>(_mu_name)),

//     _v(getADMaterialProperty<Real>("v")),

//     _e(getADMaterialProperty<Real>("e")),

//     _fp_2phase(getUserObject<HEM>("fp_2phase"))
// {
// }

// void
// ADDynamicViscosityMaterial::computeQpProperties()
// {
//   _mu[_qp] = _fp_1phase.mu_from_v_e(_v[_qp], _e[_qp]);
// }
