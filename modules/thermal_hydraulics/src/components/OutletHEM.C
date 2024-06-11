//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "OutletHEM.h"
#include "FlowModelHEM.h"

registerMooseObject("ThermalHydraulicsApp", OutletHEM);

InputParameters
OutletHEM::validParams()
{
  InputParameters params = FlowBoundaryHEM::validParams();
  params.addRequiredParam<Real>("p", "Prescribed pressure [Pa]");
  params.declareControllable("p");
  // params.addRequiredCoupledVar("alpha", "The void fraction");
  params.addClassDescription(
      "Boundary condition with prescribed pressure for 1-phase flow channels.");
  return params;
}

OutletHEM::OutletHEM(const InputParameters & params) : FlowBoundaryHEM(params)
// _alpha(coupledValue("alpha"))
{
}

void
OutletHEM::check() const
{
  FlowBoundaryHEM::check();

  auto fm = dynamic_cast<const FlowModelHEM *>(_flow_model.get());
  if (fm == nullptr)
    logError("Incompatible flow model. Make sure you use this component with HEM flow channel.");
}

void
OutletHEM::addMooseObjects()
{
  ExecFlagEnum userobject_execute_on(MooseUtils::getDefaultExecFlagEnum());
  userobject_execute_on = {EXEC_INITIAL, EXEC_LINEAR, EXEC_NONLINEAR};

  // boundary flux user object
  {
    const std::string class_name = "ADBoundaryFlux3EqnGhostPressureHEM";
    InputParameters params = _factory.getValidParams(class_name);
    params.set<Real>("p") = getParam<Real>("p");
    params.set<Real>("normal") = _normal;
    params.set<UserObjectName>("fluid_properties") = _fp_name;
    params.set<UserObjectName>("numerical_flux") = _numerical_flux_name;
    params.set<ExecFlagEnum>("execute_on") = userobject_execute_on;
    // params.set<VariableName>("alpha") = _alpha;
    getTHMProblem().addUserObject(class_name, _boundary_uo_name, params);
    connectObject(params, _boundary_uo_name, "p");
  }

  // BCs
  addWeakBC3Eqn();
}
