//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADBoundaryFlux3EqnGhostPressureHEM.h"
#include "HEM.h"
#include "THMIndices3EqnHEM.h"
#include "Numerics.h"

registerMooseObject("ThermalHydraulicsApp", ADBoundaryFlux3EqnGhostPressureHEM);

InputParameters
ADBoundaryFlux3EqnGhostPressureHEM::validParams()
{
  InputParameters params = ADBoundaryFlux3EqnGhostBase::validParams();

  params.addClassDescription("Computes boundary flux from a specified pressure for the 1-D, "
                             "1-phase, variable-area Euler equations");

  params.addRequiredParam<Real>("p", "Pressure");

  // params.addRequiredParam<MaterialPropertyName>("alpha", "Specific heat of the fluid");

  params.addRequiredParam<UserObjectName>("fluid_properties",
                                          "Name of fluid properties user object");

  params.declareControllable("p");
  return params;
}

ADBoundaryFlux3EqnGhostPressureHEM::ADBoundaryFlux3EqnGhostPressureHEM(
    const InputParameters & parameters)
  : ADBoundaryFlux3EqnGhostBase(parameters),

    _p(getParam<Real>("p")),
    // _alpha(getADMaterialProperty<Real>("alpha")),
    _fp(getUserObject<HEM>("fluid_properties"))
{
}

std::vector<ADReal>
ADBoundaryFlux3EqnGhostPressureHEM::getGhostCellSolution(const std::vector<ADReal> & U) const
{
  const ADReal rhoA = U[THM3EqnHEM::CONS_VAR_RHOA];
  const ADReal rhouA = U[THM3EqnHEM::CONS_VAR_RHOUA];
  const ADReal A = U[THM3EqnHEM::CONS_VAR_AREA];
  const ADReal alpha = U[THM3EqnHEM::CONS_VAR_ALPHA];

  const ADReal rho = rhoA / A;
  const ADReal vel = rhouA / rhoA;
  const ADReal E = _fp.e_mixture_from_p_rho(_p, rho, 0) + 0.5 * vel * vel;

  std::vector<ADReal> U_ghost(THM3EqnHEM::N_CONS_VAR);
  U_ghost[THM3EqnHEM::CONS_VAR_RHOA] = rhoA;
  U_ghost[THM3EqnHEM::CONS_VAR_RHOUA] = rhouA;
  U_ghost[THM3EqnHEM::CONS_VAR_RHOEA] = rhoA * E;
  U_ghost[THM3EqnHEM::CONS_VAR_AREA] = A;
  U_ghost[THM3EqnHEM::CONS_VAR_ALPHA] = A;

  return U_ghost;
}
