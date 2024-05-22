[FluidProperties]
  [fp_2phase]
    type = StiffenedGasTwoPhaseFluidProperties
  []
  [fp]
    type = HEM
    fp_2phase = fp_2phase
    alpha = 1
  []
[]
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
  allow_renumbering = false
[]

[Materials]
  [temperature]
    type = ADConstantMaterial
    property_name = temperature
    value = 600
  []
  [pressure]
    type = ADConstantMaterial
    property_name = pressure
    value = 15.5e+6
  []
  [const_mpropsat]
    type = ADFluidPropertiesHEM3EqnMaterial
    fp = fp
    T = temperature
    p = pressure
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
 []

[Outputs]
  csv = true
  execute_on = timestep_end
[]

