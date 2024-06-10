#Liquid Phase 0 - near minimum temperature
rhoA = 2015.0852080
rhouA = 6045.2556240
rhoEA = 7481976.7420194

# #Liquid Phase 1
# rhoA = 1437.8994142
# rhouA = 4313.6982426
# rhoEA = 1921202327.1160000

# #Liquid Phase 2
# rhoA = 1291.2327990
# rhouA = 3873.6983970
# rhoEA = 1944463933.3906100

# #Liquid Phase 3 - saturation
# rhoA = 1188.7574480
# rhouA = 3566.2723440
# rhoEA = 1906537193.1813900

# #Mixture 1 - alpha = 0.1
# rhoA = 1090.2678632
# rhouA = 3270.8035896
# rhoEA = 1765708380.3351300

# #Mixture 2 - alpha = 0.2
# rhoA = 991.7782784
# rhouA = 2975.3348352
# rhoEA = 1624879567.4888700

# #Mixture 3 - alpha = 0.3
# rhoA = 893.2886936
# rhouA = 2679.8660808
# rhoEA = 1484050754.6426100

# #Mixture 4 - alpha = 0.4
# rhoA = 794.7991088
# rhouA = 2384.3973264
# rhoEA = 1343221941.7963500

# #Mixture 5 - alpha = 0.5
# rhoA = 696.3095240
# rhouA = 2088.9285720
# rhoEA = 1202393128.9500900

# #Mixture 6 - alpha = 0.6
# rhoA = 597.8199392
# rhouA = 1793.4598176
# rhoEA = 1061564316.1038300

# #Mixture 7 - alpha = 0.7
# rhoA = 499.3303544
# rhouA = 1497.9910632
# rhoEA = 920735503.2575760

# #Mixture 8 - alpha = 0.8
# rhoA = 400.8407696
# rhouA = 1202.5223088
# rhoEA = 779906690.4113170

# #Mixture 9 - alpha = 0.9
# rhoA = 302.3511848
# rhouA = 907.0535544
# rhoEA = 639077877.5650590

# #Vapor Phase 1 - saturation
# rhoA = 203.8615446
# rhouA = 611.5846338
# rhoEA = 498248907.4233950

# #Vapor Phase 2
# rhoA = 162.8378940
# rhouA = 488.5136821
# rhoEA = 421264381.5784430

# #Vapor Phase 3 - near maximum temperature
# rhoA = 53.3211615
# rhouA = 159.9634846
# rhoEA = 214000587.2844800

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
  []
[]
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
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
    p = p
    T = T
    alpha = alpha
    fp_hem = fp
  []
[]

[Problem]
  solve = false
  # register_objects_from = 'IAPWS95App'
  # library_path = '/Users/aldelc/projects/bison_INL_desktop/iapws95/lib'
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
  [a_alpha]
    type = ADElementAverageMaterialProperty
    mat_prop = alpha
  []
  [k_press]
    type = ADElementAverageMaterialProperty
    mat_prop = p
  []
  [l_temp]
    type = ADElementAverageMaterialProperty
    mat_prop = T
  []
[]

[Outputs]
  csv = true
  execute_on = timestep_end
[]

