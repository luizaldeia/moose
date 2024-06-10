[GlobalParams]
  gravity_vector = '0 0 0'

  initial_vel = 0
  initial_p = 1e5
  initial_T = 300

  closures = simple_closures
[]

[FluidProperties]
  [fp_2phase]
    type = IAPWS95TwoPhaseFluidProperties
  []
  [fp]
    type = HEM
    fp_2phase = fp_2phase
  []
[]

[Closures]
  [simple_closures]
    type = ClosuresHEMSimple
  []
[]

[Components]
  [pipe]
    type = FlowChannelHEM
    fp = fp
    position = '0 0 0'
    orientation = '1 0 0'
    A = 1.
    length = 1
    n_elems = 10
    f = 0.001
  []

  [inlet]
    type = InletMassFlowRateTemperature1Phase
    input = 'pipe:in'
    m_dot = 1
    T = 350
  []

  [outlet]
    type = Outlet1Phase
    input = 'pipe:out'
    p = 10
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

  dt = 1e-4
  dtmin = 1.e-7

  solve_type = 'PJFNK'
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-8
  nl_max_its = 10

  l_tol = 1e-8
  l_max_its = 100

  start_time = 0.0
  num_steps = 10

  [Quadrature]
    type = GAUSS
    order = SECOND
  []
[]
