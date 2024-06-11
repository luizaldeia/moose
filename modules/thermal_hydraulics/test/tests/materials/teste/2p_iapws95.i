# #Liquid Phase 1
# rhoA = 2015.0852080
# rhouA = 6045.2556240
# rhoEA = 7481976.7420194

# #Liquid Phase 2
# rhoA = 1437.8994142
# rhouA = 4313.6982426
# rhoEA = 1921202327.1160000

# #Liquid Phase 3
# rhoA = 1291.2327990
# rhouA = 3873.6983970
# rhoEA = 1944463933.3906100

# #Liquid Phase sat
# rhoA = 1188.7574480
# rhouA = 3566.2723440
# rhoEA = 1906537193.1813900

# #Mixture Phase alpha = 0.1
# rhoA = 1090.2678577
# rhouA = 3270.8035730
# rhoEA = 1765708364.6055900

# #Mixture Phase alpha = 0.2
# rhoA = 991.7782673
# rhouA = 2975.3348020
# rhoEA = 1624879536.0297900

# #Mixture Phase alpha = 0.3
# rhoA = 893.2886770
# rhouA = 2679.8660309
# rhoEA = 1484050707.4539900

# #Mixture Phase alpha = 0.4
# rhoA = 794.7990866
# rhouA = 2384.3972599
# rhoEA = 1343221878.8781900

# #Mixture Phase alpha = 0.5
# rhoA = 696.3094963
# rhouA = 2088.9284889
# rhoEA = 1202393050.3023900

# #Mixture Phase alpha = 0.6
# rhoA = 597.8199060
# rhouA = 1793.4597179
# rhoEA = 1061564221.7265900

# #Mixture Phase alpha = 0.7
# rhoA = 499.3303156
# rhouA = 1497.9909469
# rhoEA = 920735393.1507930

# #Mixture Phase alpha = 0.8
# rhoA = 400.8407253
# rhouA = 1202.5221758
# rhoEA = 779906564.5749930

# #Mixture Phase alpha = 0.9
# rhoA = 302.3511349
# rhouA = 907.0534048
# rhoEA = 639077735.9991940

# #Vapor Phase sat
# rhoA = 203.8615446
# rhouA = 611.5846338
# rhoEA = 498248907.4233950

# #Vapor Phase 1
# rhoA = 162.8378940
# rhouA = 488.5136821
# rhoEA = 421264381.5784430

# #Vapor Phase 2
# rhoA = 53.3211615
# rhouA = 159.9634846
# rhoEA = 214000587.2844800

#teste
rhoA = 1945.46
rhouA = 5836.38
rhoEA = 7.01901e+08

[Variables]
  [A]
    initial_condition = 2
  []
  [rhoA]
    initial_condition = ${rhoA}
  []
  [rhouA]
    initial_condition = ${rhouA}
  []
  [rhoEA]
    initial_condition = ${rhoEA}
  []
[]

[FluidProperties]
  [fp_2phase]
    type = IAPWS95TwoPhaseFluidProperties
  []
  [fp]
    type = HEM
    fp_2phase = fp_2phase
    temp_guess = 300
  []
[]
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
  # ny = 1
  # nz = 1
  allow_renumbering = false
[]

[Materials]
  [const_mpropsat]
    type = ADFluidPropertiesHEM3EqnMaterial
    fp = fp
    A = A
    rhoA = rhoA
    rhouA = rhouA
    rhoEA = rhoEA
  []
  [mu]
    type = ADDynamicViscosityHEMMaterial
    v = v
    e = e
    mu = mu
    alpha = alpha
    fp_hem = fp
    T = T
    p = p
  []
[]

[Problem]
  solve = false
  # register_objects_from = 'ThermalHydraulicsApp'
  # library_path = '/Users/aldelc/projects/moose/modules/thermal_hydraulics/lib'
[]

[Executioner]
  type = Steady
[]

[Postprocessors]
  [a_rho]
    type = ADElementAverageMaterialProperty
    mat_prop = rho
  []
  [b_v]
    type = ADElementAverageMaterialProperty
    mat_prop = v
  []
  [c_e]
    type = ADElementAverageMaterialProperty
    mat_prop = e
  []
  [d_h]
    type = ADElementAverageMaterialProperty
    mat_prop = h
  []
  [e_k]
    type = ADElementAverageMaterialProperty
    mat_prop = k
  []
  [f_c]
    type = ADElementAverageMaterialProperty
    mat_prop = c
  []
  [g_cp]
    type = ADElementAverageMaterialProperty
    mat_prop = cp
  []
  [h_cv]
    type = ADElementAverageMaterialProperty
    mat_prop = cv
  []
  [i_mu]
    type = ADElementAverageMaterialProperty
    mat_prop = mu
  []
  [j_p]
    type = ADElementAverageMaterialProperty
    mat_prop = p
  []
  [k_T]
    type = ADElementAverageMaterialProperty
    mat_prop = T
  []
  # [l_vel]
  #   type = ADElementAverageMaterialProperty
  #   mat_prop = vel
  # []
  [z_alpha]
    type = ADElementAverageMaterialProperty
    mat_prop = alpha
  []
[]

[Outputs]
  csv = true
  execute_on = timestep_end
[]

