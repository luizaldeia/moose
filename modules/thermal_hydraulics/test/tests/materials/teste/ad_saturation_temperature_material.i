
# #Liquid Phase
# rhoA = 1380.730951
# rhouA = 4142.192854
# rhoEA = 1893131424.3866400

# #Mixture Phase alpha = 0.1
# rhoA = 1304.771546
# rhouA = 3914.314639
# rhoEA = 1676971885.4761200

# #Mixture Phase alpha = 0.2
# rhoA = 1173.151501
# rhouA = 3519.454503
# rhoEA = 1811548811

# #Mixture Phase alpha = 0.3
# rhoA = 1041.531455
# rhouA = 3124.594366
# rhoEA = 1749647201

# #Mixture Phase alpha = 0.4
# rhoA = 909.9114099
# rhouA = 2729.73423
# rhoEA = 1652022221

# #Mixture Phase alpha = 0.5
# rhoA = 778.2913644
# rhouA = 2334.874093
# rhoEA = 1518673871

# #Mixture Phase alpha = 0.6
# rhoA = 646.6713188
# rhouA = 1940.013956
# rhoEA = 1349602151

# #Mixture Phase alpha = 0.7
# rhoA = 515.0512733
# rhouA = 1545.15382
# rhoEA = 1144807062

#Mixture Phase alpha = 0.8
rhoA = 383.4312278
rhouA = 1150.293683
rhoEA = 618501641.3609560

# #Mixture Phase alpha = 0.9
# rhoA = 251.8111822
# rhouA = 755.4335467
# rhoEA = 628046771.6

# #Vapor Phase
# rhoA = 115.5336911
# rhouA = 346.6010734
# rhoEA = 306626936.1

# #teste
# rhoA = 1254.211225
# rhouA = 3762.633676
# rhoEA = 1880566583

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
  dim = 1
  nx = 10
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
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]

 [Postprocessors]
   [a_rho_mixture]
     type = ADElementAverageMaterialProperty
     mat_prop = rho
   []
   [b_v_mixture]
     type = ADElementAverageMaterialProperty
     mat_prop = v
   []
   [c_e_mixture]
     type = ADElementAverageMaterialProperty
     mat_prop = e
   []
   [d_h_mixture]
     type = ADElementAverageMaterialProperty
     mat_prop = h
   []
   [e_k_mixture]
     type = ADElementAverageMaterialProperty
     mat_prop = k
   []
   [f_c_mixture]
     type = ADElementAverageMaterialProperty
     mat_prop = c
   []
   [g_cp_mixture]
     type = ADElementAverageMaterialProperty
     mat_prop = cp
   []
   [h_cv_mixture]
     type = ADElementAverageMaterialProperty
     mat_prop = cv
   []
   [i_alpha_mixture]
    type = ADElementAverageMaterialProperty
    mat_prop = alpha
  []
 []

[Outputs]
  csv = true
  execute_on = timestep_end
[]

