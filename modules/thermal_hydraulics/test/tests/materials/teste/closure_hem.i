#Peak LHGR [qp_line] = [W/m]
qp_line = 47.24E+3

#Channel length [L]=[m]
L = 3.588

#Interior subchannel wetted perimeter [Pw] = [m]
Pw = '5.01325'

[GlobalParams]
  gravity_vector = '0 0 0'

  initial_vel = 3
  initial_p = 0.11e6
  initial_T = 360
  initial_alpha = 1.0
  rdg_slope_reconstruction = none
  closures = simple_closures
[]

[FluidProperties]
  [fp_2phase]
    type = IAPWS95TwoPhaseFluidProperties
  []
  # [fp_2phase]
  #   type = StiffenedGasTwoPhaseFluidProperties
  # []
  [fp]
    type = HEM
    fp_2phase = fp_2phase
    temp_guess = 360
  []
[]

[Functions]
  [heat_flux]
    type = ParsedFunction
    expression = '(q/p)*cos((pi/L)*y)'
    symbol_names = 'q p L'
    symbol_values = '${qp_line} ${Pw} ${L}'
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Components]
  [pipe]
    type = FlowChannelHEM
    fp = fp
    position = '0 0 0'
    orientation = '0 1 0'
    A = 2.
    length = 1
    n_elems = 2
    f = 0.000
  []

  [ht_pipe]
    type = HeatTransferFromHeatFluxHEM
    flow_channel = 'pipe'
    q_wall = 1E+8
    P_hf = ${Pw}
    Hw = 10000
  []

  [inlet]
    type = SolidWallHEM
    input = 'pipe:in'
  []

  # [outlet]
  #   type = SolidWallHEM
  #   input = 'pipe:out'
  # []

  [outlet]
    type = OutletHEM
    input = 'pipe:out'
    p = 0.11e+6
  []
[]

[Preconditioning]
  [pc]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  dt = 1e-3
  dtmin = 1.e-5

  solve_type = 'Newton'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-5
  nl_max_its = 100
  line_search = BASIC
  l_tol = 1e-5
  l_max_its = 100

  start_time = 0.0
  num_steps = 200

  # automatic_scaling = true
  # off_diagonals_in_auto_scaling = true
  # # compute_scaling_once = true

  resid_vs_jac_scaling_param = 0.5
[]

[Postprocessors]
  [alpha]
    type = ElementAverageValue
    variable = alpha
  []
  [T]
    type = ADElementAverageMaterialProperty
    mat_prop = T
  []
  [p]
    type = ADElementAverageMaterialProperty
    mat_prop = p
  []
[]

[Outputs]
  exodus = true
[]
