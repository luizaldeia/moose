//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RhoEAFromPressureTemperatureFunctionVelocityICHEM.h"
#include "HEM.h"
#include "Function.h"

registerMooseObject("ThermalHydraulicsApp", RhoEAFromPressureTemperatureFunctionVelocityICHEM);

InputParameters
RhoEAFromPressureTemperatureFunctionVelocityICHEM::validParams()
{
  InputParameters params = InitialCondition::validParams();
  params.addRequiredParam<UserObjectName>("fp", "The name of fluid properties object to use.");
  params.addRequiredCoupledVar("p", "The pressure");
  params.addRequiredCoupledVar("T", "The temperature");
  params.addRequiredParam<FunctionName>("vel", "The velocity");
  params.addRequiredCoupledVar("A", "Cross-sectional area");
  params.addClassDescription("Set the initial condition for rho*E*A from pressure and temperature "
                             "variables and a velocity scalar function");
  return params;
}

RhoEAFromPressureTemperatureFunctionVelocityICHEM::
    RhoEAFromPressureTemperatureFunctionVelocityICHEM(const InputParameters & parameters)
  : InitialCondition(parameters),
    _fp(getUserObject<HEM>("fp")),
    _p(coupledValue("p")),
    _T(coupledValue("T")),
    _vel(getFunction("vel")),
    _area(coupledValue("A"))
{
}

Real
RhoEAFromPressureTemperatureFunctionVelocityICHEM::value(const Point & p)
{
  const ADReal vel = _vel.value(_t, p);
  ADReal Psat, Tsat, rho_sat_l, rho_sat_g, alpha, error, rho, e;
  Real it;
  rho = 1000.0;
  it = 0.0;
  Psat = _fp.p_sat(_T[_qp]);
  Tsat = _fp.T_sat(_p[_qp]);
  do
  {

    // Calculating the saturation density for each phase
    rho_sat_l = _fp.rho_liquid_from_p_T(Psat, _T[_qp]);
    rho_sat_g = _fp.rho_vapor_from_p_T(Psat, _T[_qp]);

    // Calculating the void fraction
    alpha = (rho - rho_sat_l) / (rho_sat_g - rho_sat_l);
    // std::cout << "alpha: " << alpha.value() << "\n";
    if (alpha <= 0)
    {
      alpha = 0.0;
    }
    else if (alpha >= 1.0)
    {
      alpha = 1.0;
    }
    if (_T[_qp] > Tsat)
    {
      alpha = 1.0;
    }

    if (alpha == 0)
    {

      const ADReal rho_l = _fp.rho_liquid_from_p_T(_p[_qp], _T[_qp]);
      const ADReal e_l = _fp.e_liquid_from_p_rho(_p[_qp], rho);
      error = std::abs((rho_l - rho) / rho);
      rho = rho_l;
      e = e_l;
    }
    else if (alpha == 1)
    {
      const ADReal rho_g = _fp.rho_vapor_from_p_T(_p[_qp], _T[_qp]);
      const ADReal e_g = _fp.e_vapor_from_p_rho(_p[_qp], rho);
      error = std::abs((rho_g - rho) / rho);
      rho = rho_g;
      e = e_g;
    }

    else
    {
      const ADReal rho_m = _fp.rho_mixture_from_p_T(Psat, _T[_qp], alpha);
      const ADReal e_m = _fp.e_mixture_from_p_rho(Psat, rho, alpha);
      error = std::abs((rho_m - rho) / rho);
      rho = rho_m;
      e = e_m;
    }
    if (it >= 100)
    {
      mooseError("The maximum number of iterations was reached while trying to determine the "
                 "fluid state");
      break;
    }
    it++;
  } while (error > 1e-4);

  ADReal teste = MetaPhysicL::raw_value(rho * (e + 0.5 * vel * vel)) * _area[_qp];

  // std::cout << "\n";
  // std::cout << "========================================\n";
  // std::cout << "RhoEAFromPressureTemperatureFunctionVelocityICHEM \n";
  // std::cout << "rhoEA ---------> " << teste.value() << "\n";
  // std::cout << "========================================\n";
  // std::cout << "\n";

  return MetaPhysicL::raw_value(rho * (e + 0.5 * vel * vel)) * _area[_qp];
}
