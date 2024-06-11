//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FlowBoundaryHEM.h"

/**
 * Boundary condition with prescribed pressure for 1-phase flow channels
 */
class OutletHEM : public FlowBoundaryHEM
{
public:
  OutletHEM(const InputParameters & params);

  virtual void addMooseObjects() override;

protected:
  virtual void check() const override;

  /// The temperature
  // const VariableValue & _alpha;

public:
  static InputParameters validParams();
};
