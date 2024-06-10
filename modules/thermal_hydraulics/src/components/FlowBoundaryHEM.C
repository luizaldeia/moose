//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FlowBoundaryHEM.h"
#include "FlowChannelHEM.h"
#include "FlowModelHEM.h"

InputParameters
FlowBoundaryHEM::validParams()
{
  InputParameters params = FlowBoundary::validParams();
  return params;
}

FlowBoundaryHEM::FlowBoundaryHEM(const InputParameters & params)
  : FlowBoundary(params), _boundary_uo_name(genName(name(), "boundary_uo"))
{
}

void
FlowBoundaryHEM::init()
{
  FlowBoundary::init();

  if (hasComponentByName<FlowChannelHEM>(_connected_component_name))
  {
    const FlowChannelHEM & comp =
        getTHMProblem().getComponentByName<FlowChannelHEM>(_connected_component_name);

    _numerical_flux_name = comp.getNumericalFluxUserObjectName();
  }
}

void
FlowBoundaryHEM::check() const
{
  FlowBoundary::check();

  checkComponentOfTypeExistsByName<FlowChannelHEM>(_connected_component_name);
}

void
FlowBoundaryHEM::addWeakBC3Eqn()
{
  const std::string class_name = "ADBoundaryFlux3EqnBC";
  InputParameters params = _factory.getValidParams(class_name);
  params.set<std::vector<BoundaryName>>("boundary") = getBoundaryNames();
  params.set<Real>("normal") = _normal;
  params.set<UserObjectName>("boundary_flux") = _boundary_uo_name;
  params.set<std::vector<VariableName>>("A_linear") = {FlowModel::AREA_LINEAR};
  params.set<std::vector<VariableName>>("rhoA") = {FlowModelHEM::RHOA};
  params.set<std::vector<VariableName>>("rhouA") = {FlowModelHEM::RHOUA};
  params.set<std::vector<VariableName>>("rhoEA") = {FlowModelHEM::RHOEA};
  params.set<bool>("implicit") = getTHMProblem().getImplicitTimeIntegrationFlag();

  const std::vector<NonlinearVariableName> variables{
      FlowModelHEM::RHOA, FlowModelHEM::RHOUA, FlowModelHEM::RHOEA};

  for (const auto & var : variables)
  {
    params.set<NonlinearVariableName>("variable") = var;
    getTHMProblem().addBoundaryCondition(
        class_name, genName(name(), var, "bnd_flux_3eqn_bc"), params);
  }
}
