//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "THMPressureAuxHEM.h"
#include "HEM.h"

registerMooseObject("ThermalHydraulicsApp", THMPressureAuxHEM);

InputParameters
THMPressureAuxHEM::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("rhoA", "Conserved density");
  params.addRequiredCoupledVar("rhouA", "Conserved momentum");
  params.addRequiredCoupledVar("rhoEA", "Conserved total energy");
  params.addRequiredCoupledVar("A", "Cross-sectional area");

  params.addRequiredParam<UserObjectName>("fp", "The name of the user object for fluid properties");
  params.addClassDescription("Computes pressure given the conervative variables");
  return params;
}

THMPressureAuxHEM::THMPressureAuxHEM(const InputParameters & parameters)
  : AuxKernel(parameters),
    _rhoA(adCoupledValue("rhoA")),
    _rhouA(adCoupledValue("rhouA")),
    _rhoEA(adCoupledValue("rhoEA")),
    _area(adCoupledValue("A")),
    _fp(getUserObject<HEM>("fp"))
{
}

Real
THMPressureAuxHEM::computeValue()
{
  HEM::HEMState state = _fp.fluid_state(_rhoA[_qp], _rhouA[_qp], _rhoEA[_qp], _area[_qp]);

  ADReal teste = state.p;

  std::cout << "\n";
  std::cout << "========================================\n";
  std::cout << "THMPressureAuxHEM \n";
  std::cout << "rhoA ----------> " << _rhoA[_qp].value() << "\n";
  std::cout << "rhouA ----------> " << _rhouA[_qp].value() << "\n";
  std::cout << "rhoEA ----------> " << _rhoEA[_qp].value() << "\n";
  std::cout << "A --------------> " << _area[_qp].value() << "\n";
  std::cout << "qp -------------> " << _qp << "\n";
  std::cout << "\n";
  std::cout << "p -------------> " << teste.value() << "\n";
  std::cout << "========================================\n";
  std::cout << "\n";

  return state.p.value();
}
