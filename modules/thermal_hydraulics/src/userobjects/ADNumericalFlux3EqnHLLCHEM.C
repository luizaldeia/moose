//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADNumericalFlux3EqnHLLCHEM.h"
#include "THMIndices3EqnHEM.h"
#include "Numerics.h"

registerMooseObject("ThermalHydraulicsApp", ADNumericalFlux3EqnHLLCHEM);

InputParameters
ADNumericalFlux3EqnHLLCHEM::validParams()
{
  InputParameters params = ADNumericalFlux3EqnBase::validParams();
  params += NaNInterface::validParams();
  params.addRequiredParam<UserObjectName>("fluid_properties",
                                          "Name for fluid properties user object");
  params.addClassDescription("Computes internal side flux for the 1-D, 1-phase, variable-area "
                             "Euler equations using the HLLC approximate Riemann solver.");
  return params;
}

ADNumericalFlux3EqnHLLCHEM::ADNumericalFlux3EqnHLLCHEM(const InputParameters & parameters)
  : ADNumericalFlux3EqnBase(parameters),
    NaNInterface(this),
    _fp(getUserObject<HEM>("fluid_properties"))
{
}

void
ADNumericalFlux3EqnHLLCHEM::calcFlux(const std::vector<ADReal> & U1,
                                     const std::vector<ADReal> & U2,
                                     const ADReal & nLR_dot_d,
                                     std::vector<ADReal> & FL,
                                     std::vector<ADReal> & FR) const
{

  // extract the conserved variables and area

  const ADReal rhoA1 = U1[THM3EqnHEM::CONS_VAR_RHOA];
  const ADReal rhouA1 = U1[THM3EqnHEM::CONS_VAR_RHOUA];
  const ADReal rhoEA1 = U1[THM3EqnHEM::CONS_VAR_RHOEA];
  const ADReal A1 = U1[THM3EqnHEM::CONS_VAR_AREA];
  const ADReal ALPHA1 = U1[THM3EqnHEM::CONS_VAR_ALPHA];

  std::cout << "\n";
  std::cout << "========================================\n";
  std::cout << "ADNumericalFlux3EqnHLLCHEM \n";
  std::cout << "\n";
  std::cout << "rhoA1 ---------> " << rhoA1.value() << "\n";
  std::cout << "rhouA1 --------> " << rhouA1.value() << "\n";
  std::cout << "rhoEA1 --------> " << rhoEA1.value() << "\n";
  std::cout << "A1 ------------> " << A1.value() << "\n";
  // std::cout << "========================================\n";
  std::cout << "\n";

  const ADReal rhoA2 = U2[THM3EqnHEM::CONS_VAR_RHOA];
  const ADReal rhouA2 = U2[THM3EqnHEM::CONS_VAR_RHOUA];
  const ADReal rhoEA2 = U2[THM3EqnHEM::CONS_VAR_RHOEA];
  const ADReal A2 = U2[THM3EqnHEM::CONS_VAR_AREA];
  const ADReal ALPHA2 = U2[THM3EqnHEM::CONS_VAR_ALPHA];

  std::cout << "rhoA2 ---------> " << rhoA2.value() << "\n";
  std::cout << "rhouA2 --------> " << rhouA2.value() << "\n";
  std::cout << "rhoEA2 --------> " << rhoEA2.value() << "\n";
  std::cout << "A2 ------------> " << A2.value() << "\n";
  std::cout << "\n";

  // reference transformation normal
  const ADReal & nx = nLR_dot_d;

  // compute the primitive variables

  HEM::HEMState state1 = _fp.fluid_state(rhoA1, rhouA1, rhoEA1, A1);

  const ADReal rho1 = state1.rho;
  const ADReal u1 = state1.vel;
  const ADReal rhou1 = rho1 * u1;
  const ADReal rhoE1 = rhoEA1 / A1;
  const ADReal q1 = u1 * nx;
  const ADReal v1 = state1.v;
  const ADReal E1 = rhoEA1 / rhoA1;
  const ADReal e1 = state1.e;
  const ADReal p1 = state1.p;
  const ADReal H1 = state1.H;
  const ADReal c1 = state1.c;
  const ADReal alpha1 = state1.alpha;

  std::cout << "rho1 ----------> " << rho1.value() << "\n";
  std::cout << "u1 ------------> " << u1.value() << "\n";
  std::cout << "rhou1 ---------> " << rhou1.value() << "\n";
  std::cout << "rhoE1 ---------> " << rhoE1.value() << "\n";
  std::cout << "q1 ------------> " << q1.value() << "\n";
  std::cout << "v1 ------------> " << v1.value() << "\n";
  std::cout << "E1 ------------> " << E1.value() << "\n";
  std::cout << "e1 ------------> " << e1.value() << "\n";
  std::cout << "p1 ------------> " << p1.value() << "\n";
  std::cout << "H1 ------------> " << H1.value() << "\n";
  std::cout << "c1 ------------> " << c1.value() << "\n";
  std::cout << "alpha1 --------> " << alpha1.value() << "\n";
  std::cout << "\n";

  HEM::HEMState state2 = _fp.fluid_state(rhoA2, rhouA2, rhoEA2, A2);

  const ADReal rho2 = state2.rho;
  const ADReal u2 = state2.vel;
  const ADReal rhou2 = rho2 * u2;
  const ADReal rhoE2 = rhoEA2 / A2;
  const ADReal q2 = u2 * nx;
  const ADReal v2 = state2.v;
  const ADReal E2 = rhoEA2 / rhoA2;
  const ADReal e2 = state2.e;
  const ADReal p2 = state2.p;
  const ADReal H2 = state2.H;
  const ADReal c2 = state2.c;
  const ADReal alpha2 = state2.alpha;

  std::cout << "rho2 ----------> " << rho2.value() << "\n";
  std::cout << "u2 ------------> " << u2.value() << "\n";
  std::cout << "rhou2 ---------> " << rhou2.value() << "\n";
  std::cout << "rhoE2 ---------> " << rhoE2.value() << "\n";
  std::cout << "q2 ------------> " << q2.value() << "\n";
  std::cout << "v2 ------------> " << v2.value() << "\n";
  std::cout << "E2 ------------> " << E2.value() << "\n";
  std::cout << "e2 ------------> " << e2.value() << "\n";
  std::cout << "p2 ------------> " << p2.value() << "\n";
  std::cout << "H2 ------------> " << H2.value() << "\n";
  std::cout << "c2 ------------> " << c2.value() << "\n";
  std::cout << "alpha2 --------> " << alpha2.value() << "\n";
  std::cout << "\n";

  // compute Roe-averaged variables
  const ADReal sqrt_rho1 = std::sqrt(rho1);
  const ADReal sqrt_rho2 = std::sqrt(rho2);
  const ADReal u_roe = (sqrt_rho1 * u1 + sqrt_rho2 * u2) / (sqrt_rho1 + sqrt_rho2);
  const ADReal q_roe = u_roe * nx;
  const ADReal H_roe = (sqrt_rho1 * H1 + sqrt_rho2 * H2) / (sqrt_rho1 + sqrt_rho2);
  const ADReal h_roe = H_roe - 0.5 * u_roe * u_roe;
  const ADReal rho_roe = std::sqrt(rho1 * rho2);
  const ADReal v_roe = 1.0 / rho_roe;
  // const ADReal e_roe = _fp.e_from_v_h(v_roe, h_roe);
  // const ADReal c_roe = _fp.c_from_v_e(v_roe, e_roe);

  // https://doi.org/10.1006/jcph.2000.6515 -> Reference for this definition of the sound speed for
  // the HEM. "Such choices generally do not satisfy Property 3.1 (If F(U_L) = F(U_R),then Phi =
  // F(U_L)= F(U_R)). This property is generally thought to be a crucial aspect of Roe’s
  // numerical flux. However, our experience shows that infringing upon this property by choosing
  // other definitions for c can improve the robustness and the simplicity of the scheme without
  // serious impact on the precision."

  ADReal c_roe;

  if ((alpha1 >= 1.0) or (alpha2 >= 1.0))
  {
    c_roe = std::min(c1, c2);
  }
  else
  {
    c_roe = std::max(c1, c2);
  }

  std::cout << "u_roe ---------> " << u_roe.value() << "\n";
  std::cout << "q_roe ---------> " << q_roe.value() << "\n";
  std::cout << "H_roe ---------> " << H_roe.value() << "\n";
  std::cout << "h_roe ---------> " << h_roe.value() << "\n";
  std::cout << "rho_roe -------> " << rho_roe.value() << "\n";
  std::cout << "v_roe ---------> " << v_roe.value() << "\n";
  std::cout << "c_roe ---------> " << c_roe.value() << "\n";
  std::cout << "\n";

  // compute wave speeds
  const ADReal s1 = std::min(q1 - c1, q_roe - c_roe);
  const ADReal s2 = std::max(q2 + c2, q_roe + c_roe);
  const ADReal sm = (rho2 * q2 * (s2 - q2) - rho1 * q1 * (s1 - q1) + p1 - p2) /
                    (rho2 * (s2 - q2) - rho1 * (s1 - q1));

  std::cout << "s1 ------------> " << s1.value() << "\n";
  std::cout << "s2 ---------> " << s2.value() << "\n";
  std::cout << "sm ---------> " << sm.value() << "\n";
  std::cout << "\n";

  // compute Omega_L, Omega_R
  const ADReal omeg1 = 1.0 / (s1 - sm);
  const ADReal omeg2 = 1.0 / (s2 - sm);

  // compute p^*
  const ADReal ps = rho1 * (s1 - q1) * (sm - q1) + p1;

  // compute U_L^*, U_R^*

  const ADReal rhoLs = omeg1 * (s1 - q1) * rho1;
  const ADReal rhouLs = omeg1 * ((s1 - q1) * rhou1 + (ps - p1) * nx);
  const ADReal rhoELs = omeg1 * ((s1 - q1) * rhoE1 - p1 * q1 + ps * sm);

  std::cout << "rhoLs ---------> " << rhoLs.value() << "\n";
  std::cout << "rhouLs --------> " << rhouLs.value() << "\n";
  std::cout << "rhoELs --------> " << rhoELs.value() << "\n";
  std::cout << "\n";

  const ADReal rhoRs = omeg2 * (s2 - q2) * rho2;
  const ADReal rhouRs = omeg2 * ((s2 - q2) * rhou2 + (ps - p2) * nx);
  const ADReal rhoERs = omeg2 * ((s2 - q2) * rhoE2 - p2 * q2 + ps * sm);

  std::cout << "rhoRs ---------> " << rhoRs.value() << "\n";
  std::cout << "rhouRs --------> " << rhouRs.value() << "\n";
  std::cout << "rhoERs --------> " << rhoERs.value() << "\n";
  std::cout << "\n";

  const ADReal A_flow = computeFlowArea(U1, U2);

  std::cout << "A_flow --------> " << A_flow.value() << "\n";
  std::cout << "\n";

  // compute the fluxes
  FL.resize(THM3EqnHEM::N_EQ);
  if (s1 > 0.0)
  {
    FL[THM3EqnHEM::EQ_MASS] = u1 * rho1 * A_flow;
    FL[THM3EqnHEM::EQ_MOMENTUM] = (u1 * rhou1 + p1) * A_flow;
    FL[THM3EqnHEM::EQ_ENERGY] = u1 * (rhoE1 + p1) * A_flow;

    ADReal teste1 = FL[THM3EqnHEM::EQ_MASS];
    ADReal teste2 = FL[THM3EqnHEM::EQ_MOMENTUM];
    ADReal teste3 = FL[THM3EqnHEM::EQ_ENERGY];
    std::cout << "s1>1 \n";
    std::cout << "\n";
    std::cout << "Mass ----------> " << teste1.value() << "\n";
    std::cout << "Momentum ------> " << teste2.value() << "\n";
    std::cout << "Energy --------> " << teste3.value() << "\n";
    std::cout << "\n";

    _last_region_index = 0;
  }
  else if (s1 <= 0.0 && sm > 0.0)
  {
    FL[THM3EqnHEM::EQ_MASS] = sm * nx * rhoLs * A_flow;
    FL[THM3EqnHEM::EQ_MOMENTUM] = (sm * nx * rhouLs + ps) * A_flow;
    FL[THM3EqnHEM::EQ_ENERGY] = sm * nx * (rhoELs + ps) * A_flow;

    ADReal teste1 = FL[THM3EqnHEM::EQ_MASS];
    ADReal teste2 = FL[THM3EqnHEM::EQ_MOMENTUM];
    ADReal teste3 = FL[THM3EqnHEM::EQ_ENERGY];
    std::cout << "s1 <= 0.0 && sm > 0.0\n";
    std::cout << "\n";
    std::cout << "Mass ----------> " << teste1.value() << "\n";
    std::cout << "Momentum ------> " << teste2.value() << "\n";
    std::cout << "Energy --------> " << teste3.value() << "\n";
    std::cout << "\n";

    _last_region_index = 1;
  }
  else if (sm <= 0.0 && s2 >= 0.0)
  {
    FL[THM3EqnHEM::EQ_MASS] = sm * nx * rhoRs * A_flow;
    FL[THM3EqnHEM::EQ_MOMENTUM] = (sm * nx * rhouRs + ps) * A_flow;
    FL[THM3EqnHEM::EQ_ENERGY] = sm * nx * (rhoERs + ps) * A_flow;

    ADReal teste1 = FL[THM3EqnHEM::EQ_MASS];
    ADReal teste2 = FL[THM3EqnHEM::EQ_MOMENTUM];
    ADReal teste3 = FL[THM3EqnHEM::EQ_ENERGY];
    std::cout << "sm <= 0.0 && s2 >= 0.0\n";
    std::cout << "\n";
    std::cout << "Mass ----------> " << teste1.value() << "\n";
    std::cout << "Momentum ------> " << teste2.value() << "\n";
    std::cout << "Energy --------> " << teste3.value() << "\n";
    std::cout << "\n";

    _last_region_index = 2;
  }
  else if (s2 < 0.0)
  {
    FL[THM3EqnHEM::EQ_MASS] = u2 * rho2 * A_flow;
    FL[THM3EqnHEM::EQ_MOMENTUM] = (u2 * rhou2 + p2) * A_flow;
    FL[THM3EqnHEM::EQ_ENERGY] = u2 * (rhoE2 + p2) * A_flow;

    ADReal teste1 = FL[THM3EqnHEM::EQ_MASS];
    ADReal teste2 = FL[THM3EqnHEM::EQ_MOMENTUM];
    ADReal teste3 = FL[THM3EqnHEM::EQ_ENERGY];
    std::cout << "s2 < 0.0\n";
    std::cout << "\n";
    std::cout << "Mass ----------> " << teste1.value() << "\n";
    std::cout << "Momentum ------> " << teste2.value() << "\n";
    std::cout << "Energy --------> " << teste3.value() << "\n";
    std::cout << "\n";

    _last_region_index = 3;
  }
  else
  {
    FL = {getNaN(), getNaN(), getNaN()};
    std::cout << "nan\n";
    std::cout << "\n";
  }

  FR = FL;

  const ADReal A_wall_L = A1 - A_flow;
  FL[THM3EqnHEM::EQ_MOMENTUM] += p1 * A_wall_L;

  const ADReal A_wall_R = A2 - A_flow;
  FR[THM3EqnHEM::EQ_MOMENTUM] += p2 * A_wall_R;

  ADReal teste1 = FL[THM3EqnHEM::EQ_MASS];
  ADReal teste2 = FL[THM3EqnHEM::EQ_MOMENTUM];
  ADReal teste3 = FL[THM3EqnHEM::EQ_ENERGY];

  std::cout << "Mass ----------> " << teste1.value() << "\n";
  std::cout << "Momentum ------> " << teste2.value() << "\n";
  std::cout << "Energy --------> " << teste3.value() << "\n";
  std::cout << "========================================\n";
  std::cout << "\n";
}

ADReal
ADNumericalFlux3EqnHLLCHEM::computeFlowArea(const std::vector<ADReal> & U1,
                                            const std::vector<ADReal> & U2) const
{
  return std::min(U1[THM3EqnHEM::CONS_VAR_AREA], U2[THM3EqnHEM::CONS_VAR_AREA]);
}
