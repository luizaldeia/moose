
# #Liquid Phase
# rhoA = 1380.730951
# rhouA = 4142.192854
# rhoEA = 1893137638

# #Vapor Phase
# rhoA = 115.5336911
# rhouA = 346.6010734
# rhoEA = 306626936.1

#Mixture Phase
rhoA = 1001.171773
rhouA = 3003.51532
rhoEA = 1758037797

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
  nx = 1
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

