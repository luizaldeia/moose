//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "HEM.h"

class HEM;

/**
 * Computes velocity and thermodynamic variables from solution variables for 2-phase HEM flow.
 */
class ADFluidPropertiesHEM3EqnMaterial : public Material
{
public:
  ADFluidPropertiesHEM3EqnMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  /// Cross-sectional area
  const ADVariableValue & _area;
  const ADVariableValue & _rhoA;
  const ADVariableValue & _rhouA;
  const ADVariableValue & _rhoEA;

  /// Density
  ADMaterialProperty<Real> & _rho;

  /// Specific volume
  ADMaterialProperty<Real> & _v;

  /// Mixture specific volume
  ADMaterialProperty<Real> & _v_m;

  /// Velocity
  ADMaterialProperty<Real> & _vel;

  /// Specific internal energy
  ADMaterialProperty<Real> & _e;

  /// Mixture specific internal energy
  ADMaterialProperty<Real> & _e_m;

  /// Pressure
  ADMaterialProperty<Real> & _p;

  /// Temperature
  ADMaterialProperty<Real> & _T;

  /// Specific enthalpy
  ADMaterialProperty<Real> & _h;

  /// Specific total (stagnation) enthalpy
  ADMaterialProperty<Real> & _H;

  // /// Sound speed
  ADMaterialProperty<Real> & _c;

  /// Constant-pressure specific heat
  ADMaterialProperty<Real> & _cp;

  /// Constant-volume specific heat
  ADMaterialProperty<Real> & _cv;

  /// Thermal conductivity
  ADMaterialProperty<Real> & _k;

  /// Fluid properties
  const HEM & _fp;

public:
  static InputParameters validParams();
};
