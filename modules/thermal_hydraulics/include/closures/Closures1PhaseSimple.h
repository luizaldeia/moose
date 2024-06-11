//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Closures1PhaseBase.h"

/**
 * Simple 1-phase closures
 */
class Closures1PhaseSimple : public Closures1PhaseBase
{
public:
  Closures1PhaseSimple(const InputParameters & params);

  // ==============================================================================
  //                             Single phase flow model
  // ==============================================================================

  virtual void checkFlowChannel(const FlowChannelBase & flow_channel) const override;
  virtual void checkHeatTransfer(const HeatTransferBase & heat_transfer,
                                 const FlowChannelBase & flow_channel) const override;
  virtual void addMooseObjectsFlowChannel(const FlowChannelBase & flow_channel) override;
  virtual void addMooseObjectsHeatTransfer(const HeatTransferBase & heat_transfer,
                                           const FlowChannelBase & flow_channel) override;

  // ==============================================================================
  //                          Homogeneous Equilibrium Model
  // ==============================================================================

  virtual void checkHeatTransferHEM(const HeatTransferHEMBase & heat_transfer,
                                    const FlowChannelBase & flow_channel) const override;
  virtual void addMooseObjectsFlowChannelHEM(const FlowChannelBase & flow_channel) override;
  virtual void addMooseObjectsHeatTransferHEM(const HeatTransferHEMBase & heat_transfer,
                                              const FlowChannelBase & flow_channel) override;

protected:
  // ==============================================================================
  //                             Single phase flow model
  // ==============================================================================

  /**
   * Adds material to compute wall temperature from heat flux
   *
   * @param[in] flow_channel   Flow channel component
   */
  void addWallTemperatureFromHeatFluxMaterial(const FlowChannel1Phase & flow_channel) const;

  // ==============================================================================
  //                          Homogeneous Equilibrium Model
  // ==============================================================================

  /**
   * Adds material to compute wall temperature from heat flux
   *
   * @param[in] flow_channel   Flow channel component
   */
  void addWallTemperatureFromHeatFluxMaterialHEM(const FlowChannelHEM & flow_channel) const;

public:
  static InputParameters validParams();
};
