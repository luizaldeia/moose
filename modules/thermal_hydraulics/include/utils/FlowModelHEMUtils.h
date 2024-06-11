//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "HEM.h"
#include "THMIndices3EqnHEM.h"
#include "MooseVariable.h"

#include "libmesh/elem.h"

namespace FlowModelHEMUtils
{

/**
 * Computes the primitive solution vector from the conservative solution vector
 *
 * @param[in] U   Conservative solution vector
 * @param[in] fp  Fluid properties object
 */
template <bool is_ad>
std::vector<GenericReal<is_ad>>
computePrimitiveSolutionVector(const std::vector<GenericReal<is_ad>> & U, const HEM & fp)
{
  const auto & rhoA = U[THM3EqnHEM::CONS_VAR_RHOA];
  const auto & rhouA = U[THM3EqnHEM::CONS_VAR_RHOUA];
  const auto & rhoEA = U[THM3EqnHEM::CONS_VAR_RHOEA];
  const auto & A = U[THM3EqnHEM::CONS_VAR_AREA];
  // const auto & ALPHA = U[THM3EqnHEM::CONS_VAR_ALPHA];

  HEM::HEMState state = fp.fluid_state(rhoA, rhouA, rhoEA, A);

  const auto rho = state.rho;
  const auto vel = state.vel;
  const auto v = state.v;
  const auto e = state.e;
  const auto p = state.p;
  const auto T = state.T;
  const auto alpha = state.alpha;

  std::vector<GenericReal<is_ad>> W(THM3EqnHEM::N_PRIM_VAR);
  W[THM3EqnHEM::PRIM_VAR_PRESSURE] = p;
  W[THM3EqnHEM::PRIM_VAR_VELOCITY] = vel;
  W[THM3EqnHEM::PRIM_VAR_TEMPERATURE] = T;
  W[THM3EqnHEM::PRIM_VAR_ALPHA] = alpha;

  return W;
}

/**
 * Computes the conservative solution vector from the primitive solution vector
 *
 * @param[in] W   Primitive solution vector
 * @param[in] A   Cross-sectional area
 * @param[in] fp  Fluid properties object
 */
template <bool is_ad>
std::vector<GenericReal<is_ad>>
computeConservativeSolutionVector(const std::vector<GenericReal<is_ad>> & W,
                                  const GenericReal<is_ad> & A,
                                  const HEM & fp)
{
  const auto & p = W[THM3EqnHEM::PRIM_VAR_PRESSURE];
  const auto & T = W[THM3EqnHEM::PRIM_VAR_TEMPERATURE];
  const auto & vel = W[THM3EqnHEM::PRIM_VAR_VELOCITY];
  const auto & alpha = W[THM3EqnHEM::PRIM_VAR_ALPHA];

  ADReal rho, e, E;

  if (alpha <= 0.0)
  {
    rho = fp.rho_liquid_from_p_T(p, T);
    e = fp.e_liquid_from_p_T(p, T);
    E = e + 0.5 * vel * vel;
  }
  if (alpha >= 1.0)
  {
    rho = fp.rho_vapor_from_p_T(p, T);
    e = fp.e_vapor_from_p_T(p, T);
    E = e + 0.5 * vel * vel;
  }
  else
  {
    rho = fp.rho_mixture_from_p_T(p, T, alpha);
    e = fp.e_mixture_from_p_T(p, T, alpha);
    E = e + 0.5 * vel * vel;
  }

  std::vector<GenericReal<is_ad>> U(THM3EqnHEM::N_CONS_VAR);
  U[THM3EqnHEM::CONS_VAR_RHOA] = rho * A;
  U[THM3EqnHEM::CONS_VAR_RHOUA] = U[THM3EqnHEM::CONS_VAR_RHOA] * vel;
  U[THM3EqnHEM::CONS_VAR_RHOEA] = U[THM3EqnHEM::CONS_VAR_RHOA] * E;
  U[THM3EqnHEM::CONS_VAR_AREA] = A;

  return U;
}

/**
 * Gets the elemental conservative solution vector
 *
 * @param[in] elem         Element
 * @param[in] U_vars       Vector of conservative variable pointers
 * @param[in] is_implicit  Is implicit?
 */
template <bool is_ad>
std::vector<GenericReal<is_ad>>
getElementalSolutionVector(const Elem * elem,
                           const std::vector<MooseVariable *> & U_vars,
                           bool is_implicit)
{
  mooseAssert(elem, "The supplied element is a nullptr.");

  std::vector<GenericReal<is_ad>> U(THM3EqnHEM::N_CONS_VAR, 0.0);

  if (is_implicit)
  {
    for (unsigned int i = 0; i < THM3EqnHEM::N_CONS_VAR; i++)
    {
      mooseAssert(U_vars[i], "The supplied variable is a nullptr.");
      U[i] = U_vars[i]->getElementalValue(elem);
    }

    std::vector<dof_id_type> dof_indices;

    const std::vector<unsigned int> ind = {
        THM3EqnHEM::CONS_VAR_RHOA, THM3EqnHEM::CONS_VAR_RHOUA, THM3EqnHEM::CONS_VAR_RHOEA};
    for (unsigned int j = 0; j < ind.size(); j++)
    {
      const auto i = ind[j];
      U_vars[i]->dofMap().dof_indices(elem, dof_indices, U_vars[i]->number());
      Moose::derivInsert(U[i].derivatives(), dof_indices[0], 1.0);
    }
  }
  else
  {
    for (unsigned int i = 0; i < THM3EqnHEM::N_CONS_VAR; i++)
      U[i] = U_vars[i]->getElementalValueOld(elem);
  }

  return U;
}
}
