//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InitialCondition.h"

class HEM;

/**
 * Computes density from pressure and temperature variables
 */
class RhoFromPressureTemperatureICHEM : public InitialCondition
{
public:
  RhoFromPressureTemperatureICHEM(const InputParameters & parameters);

protected:
  virtual Real value(const Point & p);

  const HEM & _fp;
  /// The pressure
  const VariableValue & _p;
  /// The temperature
  const VariableValue & _T;

public:
  static InputParameters validParams();
};
