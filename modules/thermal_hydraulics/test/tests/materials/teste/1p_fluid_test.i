#Liquid Phase
rhoA = 1436.3915920
rhouA = 4309.1747760
rhoEA = 1828181920.3497200


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
  [fp]
    type = IAPWS95LiquidFluidProperties
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
    type = ADFluidProperties3EqnMaterial
    fp = fp
    A = A
    rhoA = rhoA
    rhouA = rhouA
    rhoEA = rhoEA
  []
  [mu]
    type = ADDynamicViscosityMaterial
    v = v
    e = e
    mu = mu
    fp_1phase = fp
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
 []

[Outputs]
  csv = true
  execute_on = timestep_end
[]

