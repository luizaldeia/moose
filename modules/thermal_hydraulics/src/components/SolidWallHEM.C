//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SolidWallHEM.h"
#include "FlowModelHEM.h"

registerMooseObject("ThermalHydraulicsApp", SolidWallHEM);

InputParameters
SolidWallHEM::validParams()
{
  InputParameters params = FlowBoundaryHEM::validParams();
  params.addClassDescription("Adds the boundary condition for a wall in single phase flow");
  return params;
}

SolidWallHEM::SolidWallHEM(const InputParameters & params) : FlowBoundaryHEM(params) {}

void
SolidWallHEM::addMooseObjects()
{
  ExecFlagEnum userobject_execute_on(MooseUtils::getDefaultExecFlagEnum());
  userobject_execute_on = {EXEC_INITIAL, EXEC_LINEAR, EXEC_NONLINEAR};

  // boundary flux user object
  {
    const std::string class_name = "ADBoundaryFlux3EqnGhostWall";
    InputParameters params = _factory.getValidParams(class_name);
    params.set<UserObjectName>("numerical_flux") = _numerical_flux_name;
    params.set<Real>("normal") = _normal;
    params.set<ExecFlagEnum>("execute_on") = userobject_execute_on;
    getTHMProblem().addUserObject(class_name, _boundary_uo_name, params);
  }

  // BCs
  addWeakBC3Eqn();
}
