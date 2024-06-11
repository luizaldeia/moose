//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatTransferHEMBase.h"
#include "FlowModelHEM.h"
#include "FlowChannelHEM.h"
#include "ClosuresBase.h"
#include "MooseUtils.h"

InputParameters
HeatTransferHEMBase::validParams()
{
  std::cout << "teste 7\n";
  InputParameters params = HeatTransferBase::validParams();
  params.addParam<FunctionName>("Hw", "Convective heat transfer coefficient [W/(m^2-K)]");
  params.declareControllable("Hw");
  return params;
}

HeatTransferHEMBase::HeatTransferHEMBase(const InputParameters & parameters)
  : HeatTransferBase(parameters)
{
  std::cout << "teste 8\n";
}

void
HeatTransferHEMBase::init()
{
  HeatTransferBase::init();
  std::cout << "teste 9\n";
}

void
HeatTransferHEMBase::initSecondary()
{
  std::cout << "teste 10\n";
  HeatTransferBase::initSecondary();

  // determine names of heat transfer variables
  if (hasComponentByName<FlowChannelHEM>(_flow_channel_name))
  {
    const FlowChannelHEM & flow_channel = getComponentByName<FlowChannelHEM>(_flow_channel_name);
    const std::string Hw_suffix = flow_channel.getHeatTransferNamesSuffix(name());

    _Hw_1phase_name = FlowModelHEM::HEAT_TRANSFER_COEFFICIENT_WALL + Hw_suffix;
    std::cout << "teste 11\n";
  }
  else
    logError("Coupled component '", _flow_channel_name, "' must be a single phase flow channel.");
}

void
HeatTransferHEMBase::check() const
{
  HeatTransferBase::check();

  if (_closures != nullptr && hasComponentByName<FlowChannelHEM>(_flow_channel_name))
    _closures->checkHeatTransfer(*this, getComponentByName<FlowChannelHEM>(_flow_channel_name));
  std::cout << "teste 12\n";
}

void
HeatTransferHEMBase::addMooseObjects()
{
  std::cout << "teste 13\n";
  HeatTransferBase::addMooseObjects();
  std::cout << "teste 13a\n";

  const auto & fc = getComponentByName<Component>(_flow_channel_name);
  std::cout << "type = " << fc.type() << std::endl;

  _closures->addMooseObjectsHeatTransfer(*this,
                                         getComponentByName<FlowChannelHEM>(_flow_channel_name));
  std::cout << "teste 13b\n";
}

const MaterialPropertyName &
HeatTransferHEMBase::getWallHeatTransferCoefficient1PhaseName() const
{
  std::cout << "teste 14\n";
  checkSetupStatus(INITIALIZED_SECONDARY);

  return _Hw_1phase_name;
}
