//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "HeatTransferFromTemperatureHEM.h"

/**
 * Heat transfer connection from a fixed temperature function for 1-phase flow
 */
class HeatTransferFromSpecifiedTemperatureHEM : public HeatTransferFromTemperatureHEM
{
public:
  HeatTransferFromSpecifiedTemperatureHEM(const InputParameters & parameters);

  virtual void addVariables() override;
  virtual void addMooseObjects() override;

protected:
  /// wall temperature function name
  const FunctionName _T_wall_fn_name;

public:
  static InputParameters validParams();
};
