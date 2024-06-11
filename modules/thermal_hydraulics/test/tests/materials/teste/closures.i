#Peak LHGR [qp_line] = [W/m]
qp_line = 47.24E+3

#Channel length [L]=[m]
L = 3.588

#Interior subchannel wetted perimeter [Pw] = [m]
Pw = '5.01325'

[GlobalParams]
  gravity_vector = '0 0 0'

  initial_vel = 3
  initial_p = 1.1e6
  initial_T = 300
  rdg_slope_reconstruction = none
  closures = simple_closures
[]

[FluidProperties]
  [fp]
    type = IAPWS95LiquidFluidProperties
  []
  # [fp_2phase]
  #   type = StiffenedGasTwoPhaseFluidProperties
  # []
[]

[Functions]
  [heat_flux]
    type = ParsedFunction
    expression = '(q/p)*cos((pi/L)*y)'
    symbol_names = 'q p L'
    symbol_values = '${qp_line} ${Pw} ${L}'
  []
  [HeatFunction]
    type = ParsedFunction
    expression = 1313127093000.32191
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Components]
  [pipe]
    type = FlowChannel1Phase
    fp = fp
    position = '0 0 0'
    orientation = '0 1 0'
    A = 2.
    length = 1
    n_elems = 15
    f = 0.001
  []

  [ht_pipe]
    type = HeatTransferFromHeatFlux1Phase
    flow_channel = 'pipe'
    q_wall = 1E+10
    P_hf = ${Pw}
    Hw = 10000
  []
  # [ht_pipe]
  #   type = HeatTransferFromSpecifiedTemperatureHEM
  #   flow_channel = pipe
  #   T_wall = 550
  #   Hw = 1.0e5
  #   P_hf = 4.4925e-2
  # []
  [inlet]
    type = SolidWall1Phase
    input = 'pipe:in'
  []

  [outlet]
    type = SolidWall1Phase
    input = 'pipe:out'
  []

  # [total_power]
  #   type = TotalPower
  #   power = 3.0e4
  # []

  # [solid]
  #   type = HeatStructureCylindrical
  #   position = '0 0 0'
  #   orientation = '0 1 0'
  #   length = 1
  #   n_elems = 10

  #   initial_T = 300

  #   names = 'fuel'
  #   widths = '0.003015'
  #   n_part_elems = '5'
  #   solid_properties = 'fuel-mat'
  #   solid_properties_T_ref = '300'
  # []
  # [hgen]
  #   type = HeatSourceFromPowerDensity
  #   hs = solid
  #   regions = 'fuel'
  #   power_density = power_density
  # []

  # [inlet]
  #   type = InletMassFlowRateTemperatureHEM
  #   input = 'pipe:in'
  #   m_dot = 1
  #   T = 500
  # []

  # [outlet]
  #   type = OutletHEM
  #   input = 'pipe:out'
  #   p = 1.1e+7
  # []
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

  dt = 1e-4
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

# [Executioner]
#   type = Transient
#   solve_type = 'newton'
#   scheme = 'bdf2'
#   line_search = 'basic'
#   abort_on_solve_fail = false

#   petsc_options = '-snes_converged_reason -ksp_converged_reason'
#   petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -mat_mffd_err
#                          -pc_factor_shift_type -pc_factor_shift_amount'
#   petsc_options_value = 'lu    mumps   1e-6          NONZERO       1e-13'
#   snesmf_reuse_base = false

#   dt = 1e-5
#   dtmin = 1e-6

#   nl_rel_tol = 1e-7
#   nl_abs_tol = 1e-4
#   nl_max_its = 50

#   l_tol = 1e-3
#   l_max_its = 300

#   automatic_scaling = true
#   off_diagonals_in_auto_scaling = true
#   compute_scaling_once = true

#   resid_vs_jac_scaling_param = 0.5

#   [Quadrature]
#     type = GAUSS
#     order = SECOND
#   []
# []
[Postprocessors]
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
