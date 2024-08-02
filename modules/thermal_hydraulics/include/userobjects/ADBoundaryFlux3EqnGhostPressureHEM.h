//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADBoundaryFlux3EqnGhostBase.h"

class HEM;

/**
 * Computes boundary flux from a specified pressure for the 1-D, 1-phase, variable-area Euler
 * equations
 */
class ADBoundaryFlux3EqnGhostPressureHEM : public ADBoundaryFlux3EqnGhostBase
{
public:
  ADBoundaryFlux3EqnGhostPressureHEM(const InputParameters & parameters);

protected:
  virtual std::vector<ADReal> getGhostCellSolution(const std::vector<ADReal> & U1) const override;

  /// Specified pressure
  const Real & _p;

  /// The void Fraction
  // const ADMaterialProperty<Real> & _alpha;

  /// Fluid properties object
  const HEM & _fp;

public:
  static InputParameters validParams();
};
