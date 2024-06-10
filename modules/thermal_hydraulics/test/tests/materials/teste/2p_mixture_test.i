# #Liquid Phase
# rhoA = 1436.3915920
# rhouA = 4309.1747760
# rhoEA = 1828181920.3497200

# #Mixture Phase alpha = 0.1
# rhoA = 1304.7715465
# rhouA = 3914.3146394
# rhoEA = 1676971885.4761200

# #Mixture Phase alpha = 0.2
# rhoA = 1173.1515009
# rhouA = 3519.4545028
# rhoEA = 1525761850.6025300

# #Mixture Phase alpha = 0.3
# rhoA = 1041.5314554
# rhouA = 3124.5943662
# rhoEA = 1374551815.7289300

# #Mixture Phase alpha = 0.4
# rhoA = 909.9114099
# rhouA = 2729.7342296
# rhoEA = 1223341780.8553400

# #Mixture Phase alpha = 0.5
# rhoA = 778.2913644
# rhouA = 2334.8740931
# rhoEA = 1072131745.9817400

# #Mixture Phase alpha = 0.6
# rhoA = 646.6713188
# rhouA = 1940.0139565
# rhoEA = 920921711.1081460

# #Mixture Phase alpha = 0.7
# rhoA = 515.0512733
# rhouA = 1545.1538199
# rhoEA = 769711676.2345510

# #Mixture Phase alpha = 0.8
# rhoA = 383.4312278
# rhouA = 1150.2936833
# rhoEA = 618501641.3609560

# #Mixture Phase alpha = 0.9
# rhoA = 251.8111822
# rhouA = 755.4335467
# rhoEA = 467291606.4873610

# #Vapor Phase
# rhoA = 120.1911367
# rhouA = 360.5734101
# rhoEA = 316081571.6137660

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
    type = StiffenedGasTwoPhaseFluidProperties
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
    alpha = alpha
    fp_2phase = fp
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
   [j_alpha]
    type = ADElementAverageMaterialProperty
    mat_prop = alpha
   []
 []

[Outputs]
  csv = true
  execute_on = timestep_end
[]

