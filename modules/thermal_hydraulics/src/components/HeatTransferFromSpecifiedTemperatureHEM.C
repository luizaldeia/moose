//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatTransferFromSpecifiedTemperatureHEM.h"

registerMooseObject("ThermalHydraulicsApp", HeatTransferFromSpecifiedTemperatureHEM);

InputParameters
HeatTransferFromSpecifiedTemperatureHEM::validParams()
{
  InputParameters params = HeatTransferFromTemperatureHEM::validParams();
  params.addRequiredParam<FunctionName>("T_wall", "Specified wall temperature [K]");
  params.declareControllable("T_wall");
  params.addClassDescription(
      "Heat transfer connection from a fixed temperature function for 1-phase flow");
  return params;
}

HeatTransferFromSpecifiedTemperatureHEM::HeatTransferFromSpecifiedTemperatureHEM(
    const InputParameters & parameters)
  : HeatTransferFromTemperatureHEM(parameters), _T_wall_fn_name(getParam<FunctionName>("T_wall"))
{
}

void
HeatTransferFromSpecifiedTemperatureHEM::addVariables()
{
  HeatTransferFromTemperatureHEM::addVariables();

  if (!_app.isRestarting())
    getTHMProblem().addFunctionIC(_T_wall_name, _T_wall_fn_name, _flow_channel_subdomains);

  makeFunctionControllableIfConstant(_T_wall_fn_name, "T_wall");
}

void
HeatTransferFromSpecifiedTemperatureHEM::addMooseObjects()
{
  HeatTransferFromTemperatureHEM::addMooseObjects();

  {
    const std::string class_name = "FunctionAux";
    InputParameters params = _factory.getValidParams(class_name);
    params.set<AuxVariableName>("variable") = _T_wall_name;
    params.set<std::vector<SubdomainName>>("block") = _flow_channel_subdomains;
    params.set<FunctionName>("function") = _T_wall_fn_name;

    ExecFlagEnum execute_on(MooseUtils::getDefaultExecFlagEnum());
    execute_on = {EXEC_INITIAL, EXEC_LINEAR};
    params.set<ExecFlagEnum>("execute_on") = execute_on;

    getTHMProblem().addAuxKernel(class_name, genName(name(), "T_wall_auxkernel"), params);
  }

  addHeatTransferKernels();
}
