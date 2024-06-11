[FluidPropertiesInterrogator]
  fp = fp
  # p = 15.5e+6
  p = 0.11e+6
  # T = 273.15
  # T = 646.9
[]

[FluidProperties]
  [fp]
    # type = IAPWS95TwoPhaseFluidProperties
    type = StiffenedGasTwoPhaseFluidProperties
    # type = IAPWS95LiquidFluidProperties
  []
[]


