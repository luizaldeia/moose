//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

class HEM;

/**
 * Computes pressure from specific volume and specific internal energy
 */
class THMPressureAuxHEM : public AuxKernel
{
public:
  static InputParameters validParams();

  THMPressureAuxHEM(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  /// density
  const ADVariableValue & _rhoA;

  /// momentum
  const ADVariableValue & _rhouA;

  /// total energy
  const ADVariableValue & _rhoEA;

  /// area
  const ADVariableValue & _area;

  const HEM & _fp;
};
