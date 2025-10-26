/*
 * Copyright (C) 2025 Giulio Barni, Eric Madge
 * This file is part of inverse_pt.
 *
 * inverse_pt is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * inverse_pt is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * inverse_pt. If not, see <https://www.gnu.org/licenses/>. 
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "inverse_pt_calculator_bindings.hpp"
#include "python_helpers.hpp"
#include "inverse_pt/profile_calculator.hpp"
#include "inverse_pt/settings.hpp"

namespace py = pybind11;
using namespace inverse_pt;
using namespace inverse_pt::profile_calculator;

void bind_profile_calculator(pybind11::module_& m) {

    // Hydrodynamic solution types
    auto hydrosol = py::enum_<HydroSolutionType>(
        m, 
        "HydroSolutionType", 
        py::arithmetic(), 
        "Types of the hydrodynamic solutions."
    )
        .value(
            "Deflagration", 
            HydroSolutionType::Deflagration,
            "weak (subsonic) deflagration"
        )
        .value(
            "Detonation",
            HydroSolutionType::Detonation,
            "weak (supersonic) detonation"
        )
        .value(
            "Hybrid",
            HydroSolutionType::Hybrid,
            "hybrid solution (supersonic deflagration)"
        )
        .value(
            "InverseDeflagration",
            HydroSolutionType::InverseDeflagration,
            "weak (supersonic) inverse deflagration"
        )
        .value(
            "InverseDetonation",
            HydroSolutionType::InverseDetonation,
            "weak (subsonic) inverse detonation"
        )
        .value(
            "InverseHybrid",
            HydroSolutionType::InverseHybrid,
            "inverse hybrid solution (subsonic inverse deflagration)"
        )
    .export_values();
    hydrosol.def("__str__", [](HydroSolutionType t) { return to_string(t); });
    hydrosol.def("__repr__", [](HydroSolutionType t) { 
        return "<HydroSolutionType: " + to_string(t) + ">"; 
    });

    // Profile calculation modes 
    py::enum_<ProfileSolveMode>(
        m, 
        "ProfileSolveMode", 
        py::arithmetic(), 
        "Selection which differential equation to solve when calculating the"
        " velocity profile."
    )
        .value("D_XI_D_V",  
            ProfileSolveMode::D_XI_D_V,  
            "always solve d(xi)/d(v)"
        )
        .value(
            "D_V_D_XI",  
            ProfileSolveMode::D_V_D_XI,  
            "always solve d(v)/d(xi)"
        )
        .value(
            "D_V_D_XI_FOR_SHOCKS",  
            ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,  
            "solve d(v)/d(xi) when the solution is expected to end in a shock"
            " and d(xi)/d(v) otherwise"
        )
    .export_values();

    m.def(
        "gamma",
        &profile_calculator::gamma,
        py::arg("v"),
        "Compute the Lorentz boost factor gamma=1/(1-v^2).\n"
        "\nParameters\n----------\n"
        "v : float\n    Velocity v\n"
        "\nReturns\n-------\n"
        "gamma : float\n    Boost factor gamma"
    );

    m.def(
        "mu",
        &profile_calculator::mu,
        py::arg("xi"), py::arg("v"),
        "Lorentz boosted velocity.\n\n"
        "This boosts the velocity v to the frame moving with velocity xi,"
        " Eq. (2.28) of [Espinosa+ (2010)].\n"
        "    mu = (xi - v)/(1 - xi v)\n"
        "\nParameters\n----------\n"
        "xi : float\n    Relative velocity xi of the new frame to which the"
        " velocity is boosted\n"
        "v : float\n    Velocity v in the original frame\n"
        "\nReturns\n-------\n"
        "mu : float\n    Boosted velocity mu(xi,v)"
    );

    m.def(
        "check_velocity",
        &check_velocity,
        py::arg("v"),
        "Check if the velocity is between 0 and 1.\n"
        "\nParameters\n----------\n"
        "v : float\n    Velocity to check\n"
        "\nReturns\n-------\n"
        "ok : bool\n"
        "    True if the velocity lies between 0 and 1, False otherwise"
    );

    m.def(
        "shock_condition",
        &shock_condition,
        py::arg("xi"), py::arg("v"),
        "Calculate the distance to the shock front.\n\n"
        "This function returns the shock condition mu(xi_sh, v_sh) xi_sh - cs^2"
        ". It returns zero if (xi,v)) corresponds to a shock front, a negative"
        " value if the point is ahead of a shock front and a positive value if"
        " it is behind the shock front.\n"
        "\nParameters\n----------\n"
        "xi : float\n    Self-similar coordinate xi = r/t\n"
        "v : float\n    Fluid vlocity v\n"
        "\nReturns\n-------\n"
        "delta : float \n    Distance to the shock front"
    );

    m.def(
        "jouguet_condition",
        &jouguet_condition,
        py::arg("xi"), py::arg("v"),
        "Calculate the distance to the Jouguet front.\n\n"
        "This function returns the Jouguet condition mu(xi_J, v_J)^2 - cs^2."
        " It returns zero if (xi,v)) corresponds to a Jouguet front, a negative"
        " value if the point is ahead of a Jouguet front and a positive value"
        " if it is behind the Jouguet front.\n"
        "\nParameters\n----------\n"
        "xi : float\n    Self-similar coordinate xi = r/t\n"
        "v : float\n    Fluid vlocity v\n"
        "\nReturns\n-------\n"
        "delta : float \n    Distance to the Jouguet front"
    );

    m.def(
        "get_v_plus_from_v_minus",
        [](double vm, double alpha) { 
            auto [v1, v2] = get_v_plus_from_v_minus(vm, alpha);
            return py::make_tuple(v1,v2); 
        },
        py::arg("vm"), py::arg("alpha"),
        "Calculate v+ from v- and alpha+.\n\n"
        "This function calculates the fluid velocity v+ ahead of the bubble (or"
        " shock) front from the velocity v- behind the front and the transition"
        " strength alpha+ ahead of the front."
        " This is Eq. (2.21) of [Espinosa+ (2010)].\n"
        "\nParameters\n----------\n"
        "vm : float\n    Fluid velocity v- behind the front\n"
        "alpha : float\n    Transition strength alpha+ ahead of the front\n"
        "\nReturns\n-------\n"
        "vp1 : float\n    The solution on the '-' branch for the fluid velocity"
        " v+ ahead of the front\n"
        "vp2 : float\n    The solution on the '+' branch for the fluid velocity"
        " v+ ahead of the front"
    );

    m.def(
        "get_v_minus_from_v_plus",
        [](double vp, double alpha) { 
            auto [v1, v2] = get_v_minus_from_v_plus(vp, alpha);
            return py::make_tuple(v1,v2); 
        },
        py::arg("vp"), py::arg("alpha"),
        "Calculate v- from v+ and alpha+.\n\n"
        "This function calculates the fluid velocity v- behind the bubble (or"
        " shock) front from the velocity v+ and the transition strength alpha+"
        " ahead of the front."
        " This is Eq. (2.21) of [Espinosa+ (2010)] solved for v-.\n"
        "\nParameters\n----------\n"
        "vp : float\n    Fluid velocity v+ ahead of the front\n"
        "alpha : float\n    Transition strength alpha+ ahead of the front\n"
        "\nReturns\n-------\n"
        "vm1 : float\n    The solution on the '-' branch for the fluid velocity"
        " v- behind the front\n"
        "vm2 : float\n    The solution on the '+' branch for the fluid velocity"
        " v- behind the front"
    );

    m.def(
        "compute_alpha_plus",
        &compute_alpha_plus,
        py::arg("vp"), py::arg("vm"),
        "Compute the transition strength from the fluid velocities in the wall"
        " rest frame.\n\n"
        "Cf. Eq. (33) of [Barni+ (2024)]\n"
        "\nParameters\n----------\n"
        "vp : float\n    Fluid velocity v+ ahead of the wall\n"
        "vm : float\n    Fluid velocity v- behind the wall\n"
        "\nReturns\n-------\n"
        "alpha : float\n"
        "    Phase transition strength alpha+ ahead of the bubble"
    );

    m.def(
        "direct_jouguet_velocity",
        &direct_jouguet_velocity,
        py::arg("alpha"),
        "Compute the Jouguet velocity for a direct phase transition, Eq. (36)"
        " of [Barni+ (2024)].\n"
        "\nParameters\n----------\n"
        "alpha : float\n    Strength parameter alpha+\n"
        "\nReturns\n-------\n" 
        "vJ : float\n    Jouguet velocity vJ(direct)"
    );

    m.def(
        "inverse_jouguet_velocity",
        &inverse_jouguet_velocity,
        py::arg("alpha"),
        "Compute the Jouguet velocity for an inverse phase transition, solution"
        " of Eq. (39) of [Barni+ (2024)].\n"
        "\nParameters\n----------\n"
        "alpha : float\n    Strength parameter alpha+\n"
        "\nReturns\n-------\n" 
        "vJinv : float\n    Inverse Jouguet velocity vJ(inverse)"
    );

    m.def(
        "alpha_max_detonation",
        &alpha_max_detonation,
        py::arg("xi_w"), 
        py::arg("return_alpha_N")=true, 
        py::arg("step_size")=1e-4,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the maximal value of alphaN or alpha+ for detonation"
        " solutions.\n\n"
        "This functions returns the maximal value of the phase transition"
        " strength alphaN far ahead of the wall (if return_alpha_N is True) or"
        " alpha+ directly in front of the wall (if return_alpha_N is False) for"
        " direct detonations.\n\n"
        "For direct detonations alpha+ and alphaN are identical, and the"
        " transition strength must be smaller than the (direct) Jouguet"
        " strength, alpha < alphaJ(direct) = alpha+(xi_w, cs).\n"
        "\nParameters\n----------\n"
        "xi_w : float\n    Wall velocity xi_w\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False).\n"
        "step_size : float, optional (default: 1e-4)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMax : float\n    The maximal value of the transition strength for"
        " direct detonations"
    );

    m.def(
        "alpha_max_deflagration",
        &alpha_max_deflagration,
        py::arg("xi_w"), 
        py::arg("return_alpha_N")=true, 
        py::arg("step_size")=1e-4,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the maximal value of alphaN or alpha+ for deflagration"
        " solutions.\n\n"
        "This functions returns the maximal value of the phase transition"
        " strength alphaN far ahead of the wall (if return_alpha_N is True) or"
        " alpha+ directly in front of the wall (if return_alpha_N is False) for"
        " direct deflagrations.\n\n"
        "For alpha+ > 1/3, the subsonic branch of the fluid velocity v+(xi_w,"
        " alpha+) ahead of the wall becomes negative.\n"
        "\nParameters\n----------\n"
        "xi_w : float\n    Wall velocity xi_w\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False).\n"
        "step_size : float, optional (default: 1e-4)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMax : float\n    The maximal value of the transition strength for"
        " direct deflagrations"
    );

    m.def(
        "alpha_max_hybrid",
        &alpha_max_hybrid,
        py::arg("xi_w"), 
        py::arg("return_alpha_N")=true, 
        py::arg("step_size")=1e-4,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the maximal value of alphaN or alpha+ for hybrid solutions."
        "\n\n"
        "This functions returns the maximal value of the phase transition"
        " strength alphaN far ahead of the wall (if return_alpha_N is True) or"
        " alpha+ directly in front of the wall (if return_alpha_N is False) for"
        " direct hybrids.\n\n"
        "For alpha+ > 1/3, the subsonic branch of the fluid velocity v+(cs,"
        " alpha+) ahead of the wall becomes negative.\n"
        "\nParameters\n----------\n"
        "xi_w : float\n    Wall velocity xi_w\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False).\n"
        "step_size : float, optional (default: 1e-4)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMax : float\n    The maximal value of the transition strength for"
        " direct hybrids"
    );

    m.def(
        "alpha_min_hybrid",
        &alpha_min_hybrid,
        py::arg("xi_w"), 
        py::arg("return_alpha_N")=true, 
        py::arg("step_size")=1e-4,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the minimal value of alphaN or alpha+ for hybrid solutions."
        "\n\n"
        "This functions returns the minimal value of the phase transition"
        " strength alphaN far ahead of the wall (if return_alpha_N is True) or"
        " alpha+ directly in front of the wall (if return_alpha_N is False) for"
        " direct hybrids.\n\n"
        "In order to be able to evolve towards a shock front, the matching"
        " solution at the bubble wall needs to satisfy xi_w v_+(cs,alpha+) <"
        " cs^2, hence alpha+ > (xi_w - cs)^2/(3 xi_w^2 - cs^2).\n"
        "\nParameters\n----------\n"
        "xi_w : float\n    Wall velocity xi_w\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False).\n"
        "step_size : float, optional (default: 1e-4)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin : float\n    The minimal value of the transition strength for"
        " direct hybrids"
    );

    m.def(
        "alpha_min_inverse_deflagration",
        &alpha_min_inverse_deflagration,
        py::arg("xi_w"), 
        py::arg("return_alpha_N")=true, 
        py::arg("step_size")=1e-4,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the minimal value of alphaN or alpha+ for inverse deflagration"
        " solutions.\n\n"
        "This functions returns the minimal value of the phase transition"
        " strength alphaN far ahead  of the wall (if return_alpha_N is True) or"
        " alpha+ directly in front of the wall (if return_alpha_N is False) for"
        " inverse deflagrations.\n\n"
        "For inverse deflagrations alpha+ and alphaN are identical."
        " For alpha+ < (cs^2 - xi_w)/(1 + xi_w), the supersonic branch of the"
        " fluid velocity v-(xi_w, alpha+) behind the wall becomes superluminal."
        "\n"
        "\nParameters\n----------\n"
        "xi_w : float\n    Wall velocity xi_w\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False).\n"
        "step_size : float, optional (default: 1e-4)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin: float\n    The minimal value of the transition strength for"
        " inverse deflagrations"
    );

    m.def(
        "alpha_min_inverse_hybrid",
        &alpha_min_inverse_hybrid,
        py::arg("xi_w"), 
        py::arg("return_alpha_N")=true, 
        py::arg("step_size")=1e-4,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the minimal value of alphaN or alpha+ for inverse hybrid"
        " solutions.\n\n"
        "This functions returns the minimal value of the phase transition"
        " strength alphaN far ahead of the wall (if return_alpha_N is True) or"
        " alpha+ directly in front of the wall (if return_alpha_N is False) for"
        " inverse hybrids.\n\n"
        " For alpha+ < 1 - 2 cs, the supersonic branch of the fluid velocity"
        " v-(cs, alpha+) behind the wall becomes superluminal.\n"
        "\nParameters\n----------\n"
        "xi_w : float\n    Wall velocity xi_w\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False).\n"
        "step_size : float, optional (default: 1e-4)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin: float\n    The minimal value of the transition strength for"
        " inverse hybrids"
    );

    m.def(
        "alpha_min_inverse_detonation",
        &alpha_min_inverse_detonation,
        py::arg("xi_w"), 
        py::arg("return_alpha_N")=true, 
        py::arg("step_size")=1e-4,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the minimal value of alphaN or alpha+ for inverse detonation"
        " solutions.\n\n"
        "This functions returns the minimal value of the phase transition"
        " strength alphaN far ahead of the wall (if return_alpha_N is True) or"
        " alpha+ directly in front of the wall (if return_alpha_N is False) for"
        " inverse detonations.\n\n"
        "For inverse detonations, the transition strength must be larger than"
        " the inverse Jouguet strength, alpha+ > alphaJ(inv) = alpha+(c_s,"
        " xi_w).\n"
        "\nNote\n----\n"
        "This is also the maximal strength for inverse hybrids.\n"
        "\nParameters\n----------\n"
        "xi_w : float\n    Wall velocity xi_w\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False).\n"
        "step_size : float, optional (default: 1e-4)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin: float\n    The minimal value of the transition strength for"
        " inverse detonations"
    );

    m.def(
        "alpha_max",
        py::overload_cast<
            HydroSolutionType, double, bool, double, ProfileSolveMode
        >(&alpha_max),
        py::arg("sol_type"), 
        py::arg("xi_w"), 
        py::arg("return_alpha_N")=true, 
        py::arg("step_size")=1e-4,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the maximal value of alphaN or alpha+ for the given solution"
        " type.\n\n"
        "This functions returns the maximal value of the phase transition"
        " strength alphaN far ahead of the wall (if return_alpha_N is True) or"
        " alpha+ directly in front of the wall (if return_alpha_N is False) for"
        " the solution type specified in sol_type.\n"
        "\nParameters\n----------\n"
        "sol_type : HydroSolutionType\n    Type of the solution\n"
        "xi_w : float\n    Wall velocity xi_w\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False).\n"
        "step_size : float, optional (default: 1e-4)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMax : float\n    The maximal value of the transition strength for"
        " the specified solution type"
    );

    m.def(
        "alpha_min",
        py::overload_cast<
            HydroSolutionType, double, bool, double, ProfileSolveMode
        >(&alpha_min),
        py::arg("sol_type"), 
        py::arg("xi_w"), 
        py::arg("return_alpha_N")=true, 
        py::arg("step_size")=1e-4,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the minimal value of alphaN or alpha+ for the given solution"
        " type.\n\n"
        "This functions returns the minimal value of the phase transition"
        " strength alphaN far ahead of the wall (if return_alpha_N is True) or"
        " alpha+ directly in front of the wall (if return_alpha_N is False) for"
        " the solution type specified in sol_type.\n"
        "\nParameters\n----------\n"
        "sol_type : HydroSolutionType\n    Type of the solution\n"
        "xi_w : float\n    Wall velocity xi_w\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False).\n"
        "step_size : float, optional (default: 1e-4)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin : float\n    The minimal value of the transition strength for"
        " the specified solution type"
    );

    m.def(
        "alpha_max",
        py::overload_cast<double, bool, double, ProfileSolveMode>(&alpha_max),
        py::arg("xi_w"), 
        py::arg("return_alpha_N")=true, 
        py::arg("step_size")=1e-4,
        py::arg("solve_mode") = ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the maximal value of alphaN or alpha+.\n\n"
        "This functions returns the maximal value of the phase transition"
        " strength alphaN far ahead of the wall (if return_alpha_N is True) or"
        " alpha+ directly in front of the wall (if return_alpha_N is False).\n"
        "\nParameters\n----------\n"
        "xi_w : float\n    Wall velocity xi_w\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False).\n"
        "step_size : float, optional (default: 1e-4)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMax : float\n    The maximal value of the transition strength"
    );

    m.def(
        "alpha_min",
        py::overload_cast<double, bool, double, ProfileSolveMode>(&alpha_min),
        py::arg("xi_w"), 
        py::arg("return_alpha_N")=true, 
        py::arg("step_size")=1e-4,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the minimal value of alphaN or alpha+.\n\n"
        "This functions returns the minimal value of the phase transition"
        " strength alphaN far ahead of the wall (if return_alpha_N is True) or"
        " alpha+ directly in front of the wall (if return_alpha_N is False).\n"
        "\nParameters\n----------\n"
        "xi_w : float\n    Wall velocity xi_w\n"
        "return_alpha_N : bool, optional (default: True)\n    Whether the"
        " function should return max(alphaN) (True) or max(alpha+) (False).\n"
        "step_size : float, optional (default: 1e-4)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n" 
        "alphaMin : float\n    The minimal value of the transition strength"
    );

    m.def(
        "classify_wave_fluid_frame",
        &classify_wave_fluid_frame,
        py::arg("xi_wall"), 
        py::arg("alpha"), 
        py::arg("step_size")=1e-4,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Classify the type of the solution based on v_w and alphaN in the fluid"
        " frame.\n"
        "\nParameters\n----------\n"
        "xi_wall : float\n    Bubble wall velocity xi_w\n"
        "alpha : float\n"
        "    Phase transition strength alphaN (far ahead of the wall)\n"
        "step_size : float, optional (default: 1e-4)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "hydrotype : HydroSolutionType\n    Type of the hydrodynamic solution"
    );

    m.def(
        "classify_wave_wall_frame",
        &classify_wave_wall_frame,
        py::arg("vp"), py::arg("vm"),
        "Classify the type of the solution based on v+ and v- in the wall rest"
        " frame.\n"
        "\nParameters\n----------\n"
        "vp : float\n    Fluid velocity v+ ahead of the wall\n"
        "vm : float\n    Fluid velocity v- behind the wall\n"
        "\nReturns\n-------\n"
        "hydrotype : HydroSolutionType\n    Type of the hydrodynamic solution"
    );

    m.def(
        "compute_fluid_velocities",
        [](double xi_wall, double alpha, HydroSolutionType sol_type) { 
            auto [vp, vm] = compute_fluid_velocities(
                xi_wall, alpha, sol_type
            );
            return py::make_tuple(vp,vm); 
        },
        py::arg("xi_wall"), py::arg("alpha"), py::arg("sol_type"),
        "Compute the fluid velocities v+ and v- in the bubble rest frame.\n"
        "\nNote\n----\n"
        "For direct and inverse hybrids, v- (direct) and v+ (inverse) are fixed"
        " to cs, respectively."
        " This function then returns xi_w instead of cs.\n"
        "\nParameters\n----------\n"
        "xi_wall : float\n    Bubble wall velocity xi_w\n"
        "alpha : float\n"
        "    Phase transition strength alpha+ ahead of the wall\n"
        "sol_type : HydroSolutionType\n    Type of the solution.\n"
        "\nReturns\n-------\n"
        "vp : float\n"
        "    The fluid velocity v+ ahead of the wall in the fluid rest frame\n"
        "vm : float\n"
        "    The fluid velocity v- behind the wall in the fluid rest frame"
    );

    m.def(
        "dxi_dv",
        &dxi_dv,
        py::arg("v"), py::arg("xi"),
        "Derivative dxi/dv.\n\n"
        "We solve this differential equation to obtain xi(v) and then invert to"
        " obtain v(xi). This is Eq. (2.27) of [Espinosa+ (2010)].\n"
        "\nParameters\n----------\n"
        "v : float\n    Fluid velocity v\n"
        "xi : float\n    self-similar coordinate xi=r/t\n"
        "\nReturns\n-------\n"
        "dxi_dv : float\n    Derivative d(xi)/d(v)"
    );

    m.def(
        "dv_dxi",
        &dv_dxi,
        py::arg("xi"), py::arg("v"),
        "Derivative dv/dxi.\n\n"
        "We solve this differential equation to directly obtain v(xi)."
        " This is Eq. (2.27) of [Espinosa+ (2010)].\n"
        "\nParameters\n----------\n"
        "xi : float\n    self-similar coordinate xi=r/t\n"
        "v : float\n    Fluid velocity v\n"
        "\nReturns\n-------\n"
        "dv_dxi : float\n    Derivative d(v)/d(xi)"
    );

    m.def(
        "integrate_xi_of_v",
        [](
            py::array_t<double> xi_arr, 
            py::array_t<double> v_arr, 
            double v_stop, 
            double dv, 
            bool stop_at_shock
        ) {
            auto xi_vec = np2vec(xi_arr);
            auto v_vec = np2vec(v_arr);
            int status = integrate_xi_of_v(
                xi_vec, v_vec, v_stop, dv, stop_at_shock
            );
            return py::make_tuple(status, vec2np(xi_vec), vec2np(v_vec));
        },
        py::arg("xi_arr"), 
        py::arg("v_arr"), 
        py::arg("v_stop"), 
        py::arg("dv"), 
        py::arg("stop_at_shock")=true,
        "Solve the differential equation d(xi)/d(v) for the velocity profile."
        "\n\n"
        "This function solves the differential equation using a fourth-order"
        " Runge-Kutta method. The first two arguments are arrays to which the"
        " solutions are appended. These arrays must already contain the initial"
        " condition.\n"
        "\nParameters\n----------\n"
        "xi_arr : numpy.ndarray, shape (N,)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v_arr : numpy.ndarray, shape (N,)\n"
        "    Array containing the fluid velocity v\n"
        "v_stop : float\n    Final fluid velocity vf\n"
        "dv : float\n    Fluid velocity step size\n"
        "stop_at_shock : bool, optional (default: True)\n"
        "    Whether to stop the evolution when a shock front is encoutered\n"
        "\nReturns\n-------\n"
        "status : int\n"
        "    A status code indicating the reason for termination:\n"
        "        0  - reached v_stop\n"
        "        1  - encountered shock\n"
        "        -1 - entered unphysical region for xi\n"
        "        -2 - reached superluminal fluid velocity v\n"
        "xi_arr : numpy.ndarray, shape (N',)\n  "
        "  Array with the self-similar coordinate xi of the solution appended\n"
        "v_arr : numpy.ndarray, shape (N',)\n"
        "    Array with the fluid velocity v of the solution appended"
    );

    m.def(
        "integrate_v_of_xi",
        [](
            py::array_t<double> xi_arr, 
            py::array_t<double> v_arr, 
            double v_stop, 
            double dxi, 
            bool stop_at_shock
        ) {
            auto xi_vec = np2vec(xi_arr);
            auto v_vec = np2vec(v_arr);
            int status = integrate_v_of_xi(
                xi_vec, v_vec, v_stop, dxi, stop_at_shock
            );
            return py::make_tuple(status, vec2np(xi_vec), vec2np(v_vec));
        },
        py::arg("xi_arr"), 
        py::arg("v_arr"), 
        py::arg("v_stop"), 
        py::arg("dxi"), 
        py::arg("stop_at_shock")=true,
        "Solve the differential equation d(v)/d(xi) for the velocity profile."
        "\n\n"
        "This function solves the differential equation using a fourth-order"
        " Runge-Kutta method. The first two arguments are arrays to which the"
        " solutions are appended. These arrays must already contain the initial"
        " condition.\n"
        "\nParameters\n----------\n"
        "xi_arr : numpy.ndarray, shape (N,)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v_arr : numpy.ndarray, shape (N,)\n"
        "    Array containing the fluid velocity v\n"
        "v_stop : float\n    Final fluid velocity vf\n"
        "dxi : float\n    Self-similar coordinate step size\n"
        "stop_at_shock : bool, optional (default: True)\n"
        "    Whether to stop the evolution when a shock front is encoutered\n"
        "\nReturns\n-------\n"
        "status : int\n"
        "    A status code indicating the reason for termination:\n"
        "        0  - reached v_stop\n"
        "        1  - encountered shock\n"
        "        -1 - entered unphysical region for xi\n"
        "        -2 - reached superluminal fluid velocity v\n"
        "        -3 - crossed Jouguet front\n"
        "xi_arr : numpy.ndarray, shape (N',)\n  "
        "  Array with the self-similar coordinate xi of the solution appended\n"
        "v_arr : numpy.ndarray, shape (N',)\n"
        "    Array with the fluid velocity v of the solution appended"
    );

    m.def(
        "integrate_velocity_profile",
        [](
            py::array_t<double> xi_arr, 
            py::array_t<double> v_arr, 
            double v_stop, 
            double step_size, 
            bool stop_at_shock,
            ProfileSolveMode solve_mode
        ) {
            auto xi_vec = np2vec(xi_arr);
            auto v_vec = np2vec(v_arr);
            int status = integrate_velocity_profile(
                xi_vec, v_vec, v_stop, step_size, stop_at_shock, solve_mode
            );
            return py::make_tuple(status, vec2np(xi_vec), vec2np(v_vec));
        },
        py::arg("xi_arr"), 
        py::arg("v_arr"), 
        py::arg("v_stop"), 
        py::arg("step_size"), 
        py::arg("stop_at_shock")=true,
        py::arg("solve_mode")=ProfileSolveMode::D_XI_D_V,
        "Solve the differential equation for the velocity profile.\n\n"
        "This function solves the differential equation using a fourth-order"
        " Runge-Kutta method. The first two arguments are arrays to which the"
        " solutions are appended. These arrays must already contain the initial"
        " condition.\n"
        "\nParameters\n----------\n"
        "xi_arr : numpy.ndarray, shape (N,)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v_arr : numpy.ndarray, shape (N,)\n"
        "    Array containing the fluid velocity v\n"
        "v_stop : float\n    Final fluid velocity vf\n"
        "step_size : float\n    Fluid velocity or self-similar coordinate step"
        " size (depending on solve_dxi_dv)\n"
        "stop_at_shock : bool, optional (default: True)\n"
        "    Whether to stop the evolution when a shock front is encoutered\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_XI_D_V)\n    ProfileSolveMode switch determining which differential"
        " equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "status : int\n"
        "    A status code indicating the reason for termination:\n"
        "        0  - reached v_stop\n"
        "        1  - encountered shock\n"
        "        -1 - entered unphysical region for xi\n"
        "        -2 - reached superluminal fluid velocity v\n"
        "        -3 - crossed Jouguet front\n"
        "xi_arr : numpy.ndarray, shape (N',)\n  "
        "  Array with the self-similar coordinate xi of the solution appended\n"
        "v_arr : numpy.ndarray, shape (N',)\n"
        "    Array with the fluid velocity v of the solution appended"
    );    

    m.def(
        "compute_velocity_profile_direct_detonation",
        [](
            double v_plus, 
            double v_minus, 
            double step_size, 
            ProfileSolveMode solve_mode
        ) {
            std::vector<double> xi, v;
            size_t idx_wall = compute_velocity_profile_direct_detonation(
                xi, v, v_plus, v_minus, step_size, solve_mode
            );
            return py::make_tuple(idx_wall, vec2np(xi), vec2np(v));
        },
        py::arg("v_plus"), 
        py::arg("v_minus"), 
        py::arg("step_size")=1e-3,
        py::arg("solve_mode") = ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the velocity profile for a direct detonation.\n\n"
        "This functions solves the fluid velocity profile for a direct"
        " detonation and returns the results of the velocity v as a function of"
        " the self-similar coordiate xi.\n"
        "\nParameters\n----------\n"
        "v_plus : float\n    Fluid velocity v+ in front of the wall\n"
        "v_minus : float\n    Fluid velocity v- behind the wall\n"
        "step_size : float, optional (default: 0.001)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : integer\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "     Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "     Array containing the fluid velocity v"
    );

    m.def(
        "compute_velocity_profile_direct_deflagration",
        [](
            double v_plus, 
            double v_minus, 
            double step_size, 
            ProfileSolveMode solve_mode
        ) {
            std::vector<double> xi, v;
            size_t idx_wall = compute_velocity_profile_direct_deflagration(
                xi, v, v_plus, v_minus, step_size, solve_mode
            );
            return py::make_tuple(idx_wall, vec2np(xi), vec2np(v));
        },
        py::arg("v_plus"), 
        py::arg("v_minus"), 
        py::arg("step_size")=1e-3,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the velocity profile for a direct deflagration.\n\n"
        "This functions solves the fluid velocity profile for a direct"
        " deflagration and returns the results of the velocity v as a function"
        " of the self-similar coordiate xi.\n"
        "\nParameters\n----------\n"
        "v_plus : float\n    Fluid velocity v+ in front of the wall\n"
        "v_minus : float\n    Fluid velocity v- behind the wall\n"
        "step_size : float, optional (default: 0.001)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : integer\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "     Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "    Array containing the fluid velocity v"
    );

    m.def(
        "compute_velocity_profile_direct_hybrid",
        [](
            double v_plus, 
            double xi_wall, 
            double step_size, 
            ProfileSolveMode solve_mode
        ) {
            std::vector<double> xi, v;
            size_t idx_wall = compute_velocity_profile_direct_hybrid(
                xi, v, v_plus, xi_wall, step_size, solve_mode
            );
            return py::make_tuple(idx_wall, vec2np(xi), vec2np(v));
        },
        py::arg("v_plus"), 
        py::arg("v_minus"), 
        py::arg("step_size")=1e-3,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the velocity profile for a direct hybrid.\n\n"
        "This functions solves the fluid velocity profile for a direct hybrid"
        " and returns the results of the velocity v as a function of the"
        " self-similar coordiate xi.\n"
        "\nParameters\n----------\n"
        "v_plus : float\n    Fluid velocity v+ in front of the wall\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "step_size : float, optional (default: 0.001)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : integer\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "     Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "     Array containing the fluid velocity v"
    );

    m.def(
        "compute_velocity_profile_inverse_detonation",
        [](
            double v_plus, 
            double v_minus, 
            double step_size, 
            ProfileSolveMode solve_mode
        ) {
            std::vector<double> xi, v;
            size_t idx_wall = compute_velocity_profile_inverse_detonation(
                xi, v, v_plus, v_minus, step_size, solve_mode
            );
            return py::make_tuple(idx_wall, vec2np(xi), vec2np(v));
        },
        py::arg("v_plus"), 
        py::arg("v_minus"), 
        py::arg("step_size")=1e-3,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the velocity profile for an inverse detonation.\n\n"
        "This functions solves the fluid velocity profile for an inverse"
        " detonation and returns the results of the velocity v as a function of"
        " the self-similar coordiate xi.\n"
        "\nParameters\n----------\n"
        "v_plus : float\n    Fluid velocity v+ in front of the wall\n"
        "v_minus : float\n    Fluid velocity v- behind the wall\n"
        "step_size : float, optional (default: 0.001)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : integer\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "     Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "     Array containing the fluid velocity v"
    );

    m.def(
        "compute_velocity_profile_inverse_deflagration",
        [](
            double v_plus, 
            double v_minus, 
            double step_size, 
            ProfileSolveMode solve_mode
        ) {
            std::vector<double> xi, v;
            size_t idx_wall = compute_velocity_profile_inverse_deflagration(
                xi, v, v_plus, v_minus, step_size, solve_mode
            );
            return py::make_tuple(idx_wall, vec2np(xi), vec2np(v));
        },
        py::arg("v_plus"), 
        py::arg("v_minus"), 
        py::arg("step_size")=1e-3,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the velocity profile for an inverse deflagration.\n\n"
        "This functions solves the fluid velocity profile for an inverse"
        " deflagration and returns the results of the velocity v as a function"
        " of the self-similar coordiate xi.\n"
        "\nParameters\n----------\n"
        "v_plus : float\n    Fluid velocity v+ in front of the wall\n"
        "v_minus : float\n    Fluid velocity v- behind the wall\n"
        "step_size : float, optional (default: 0.001)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : integer\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "     Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "     Array containing the fluid velocity v"
    );

    m.def(
        "compute_velocity_profile_inverse_hybrid",
        []( 
            double xi_wall, 
            double v_minus, 
            double step_size, 
            ProfileSolveMode solve_mode
        ) {
            std::vector<double> xi, v;
            size_t idx_wall = compute_velocity_profile_inverse_hybrid(
                xi, v, xi_wall, v_minus, step_size, solve_mode
            );
            return py::make_tuple(idx_wall, vec2np(xi), vec2np(v));
        },
        py::arg("xi_wall"), 
        py::arg("v_minus"), 
        py::arg("step_size")=1e-3,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the velocity profile for an inverse hybrid.\n\n"
        "This functions solves the fluid velocity profile for an inverse hybrid"
        " and returns the results of the velocity v as a function of the"
        " self-similar coordiate xi.\n"
        "\nParameters\n----------\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "v_minus : float\n    Fluid velocity v- behind the wall\n"
        "step_size : float, optional (default: 0.001)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : integer\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "     Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "     Array containing the fluid velocity v"
    );

    m.def(
        "get_shock_index",
        // throws exception if solution has no shock front -> no idx2obj needed
        get_shock_index, 
        py::arg("sol_type"), py::arg("size"),
        "Get the position of the shock front.\n"
        "\nParameters\n----------\n"
        "sol_type : HydroSolutionType\n    Type of the solution.\n"
        "size : integer\n    Length of the profile vectors\n"
        "\nReturns\n-------\n"
        "idx : integer\n    Index of the first element outside the shock front"
    );

    m.def(
        "get_shock_front",
        [](
            profile_calculator::HydroSolutionType sol_type, 
            py::array_t<double> xi, 
            py::array_t<double> v
        ) { 
            auto [xish, vsh] = get_shock_front(sol_type, np2vec(xi), np2vec(v));
            return py::make_tuple(xish, vsh);
        },
        py::arg("sol_type"), py::arg("xi"), py::arg("v"),
        "Get the coordinate and velocity of the shock front.\n"
        "\nParameters\n----------\n"
        "sol_type : HydroSolutionType\n    Type of the solution.\n"
        "xi : numpy.ndarray, shape (N,)\n"
        "    Self-similar coordinate xi=r/t of the profile\n"
        "v : numpy.ndarray, shape (N,)\n    Fluid velocity profile v(xi)\n"
        "\nReturns\n-------\n"
        "xi_sh : float\n    self-similar coordinate of the shock front\n"
        "v_sh : float\n"
        "    fluid velocity at the shock front (in the fluid rest frame)"
    );

    m.def(
        "compute_velocity_profile",
        [](
            double xi_wall, 
            double alpha_N, 
            double step_size, 
            ProfileSolveMode solve_mode
         ) {
            std::vector<double> xi, v;
            size_t idx_wall = compute_velocity_profile(
                xi, v, xi_wall, alpha_N, step_size, solve_mode
            );
            return py::make_tuple(idx_wall, vec2np(xi), vec2np(v));
        },
        py::arg("xi_wall"), 
        py::arg("alpha_N"), 
        py::arg("step_size")=1e-3,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the velocity profile for a bubble in a cosmological phase"
        " transition.\n\n"
        "This functions solves the fluid velocity profile for a cosmological"
        " phase transition and returns the results of the velocity v as a"
        " function of the self-similar coordiate xi.\n"
        "\nParameters\n----------\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "alpha_N : float\n"
        "    Transition strength alphaN (far in front of the wall)\n"
        "step_size : float, optional (default: 0.001)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : integer\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "     Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "     Array containing the fluid velocity v"
    );

    m.def(
        "compute_velocity_profile",
        [](
            HydroSolutionType sol_type, 
            double xi_wall, 
            double alpha_p, 
            double step_size, 
            ProfileSolveMode solve_mode
        ) {
            std::vector<double> xi, v;
            size_t idx_wall = compute_velocity_profile(
                sol_type, xi, v, xi_wall, alpha_p, step_size, solve_mode
            );
            return py::make_tuple(idx_wall, vec2np(xi), vec2np(v));
        },
        py::arg("sol_type"), 
        py::arg("xi_wall"),
        py::arg("alpha_p"), 
        py::arg("step_size")=1e-3,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the velocity profile for a bubble in a cosmological phase"
        " transition.\n\n"
        "This functions solves the fluid velocity profile for a cosmological"
        " phase transition and returns the results of the velocity v as a"
        " function of the self-similar coordiate xi.\n"
        "\nParameters\n----------\n"
        "sol_type : HydroSolutionType\n    Type of the solution.\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "alpha_p : float\n    Transition strength alpha+ ahead of the wall\n"
        "step_size : float, optional (default: 0.001)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "idx_wall : integer\n    Index of the first element outside the bubble"
        " (i.e. at the wall velocity xi_w)\n"
        "xi : numpy.ndarray, shape (N',)\n"
        "     Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N',)\n"
        "     Array containing the fluid velocity v"
    );

    m.def(
        "integrate_enthalpy_segement",
        &integrate_enthalpy_segment,
        py::arg("xi1"), py::arg("xi2"), py::arg("v1"), py::arg("v2"),
        "Calculate a segment of the enthalpy density integral.\n\n"
        "This function calculates the segment of the integral in Eq. (2.29) of"
        " [Espinosa+ (2010)] assuming that v is linear in the segment (xi1,xi2)"
        ", i.e. xi = xi0 + v d(xi)/d(v), where xi0 = (xi1 v2 - xi2 v1)/"
        "(v2 - v1) and d(xi)/d(v) = (xi2-xi1)/(v2-v1). The integral can the be"
        " calculated using partial fraction decomposition.\n"
        "\nParameters\n----------\n"
        "xi1 : float\n"
        "    Self-similar coordinate xi1 at the beginning of the segment\n"
        "xi2 : float\n"
        "    Self-similar coordinate xi2 at the end of the segment\n"
        "v1 : float\n"
        "    Fluid velocity v1 at the beginning of the segment\n"
        "v2 : float\n"
        "    Fluid velocity coordinate v2 at the end of the segment\n"
        "\nReturns\n-------\n"
        "wsegment : float\n    Enthalpy density within the segment (x1,x2)."
    );

    m.def(
        "compute_enthalpy_density",
        [](
            py::array_t<double> xi, 
            py::array_t<double> v, 
            size_t idx_wall, 
            py::object idx_shock
        ) { 
            return vec2np(compute_enthalpy_density(
                np2vec(xi), np2vec(v), idx_wall, obj2idx(idx_shock)
            ));
        },
        py::arg("xi"), 
        py::arg("v"), 
        py::arg("idx_wall"), 
        py::arg("idx_shock")=py::none(),
        "Compute the enthalpy density profile.\n\n"
        "Compute the enthalpy density profile from the fluid velocity profile,"
        " Eq. (2.29) of [Espinosa+ (2010)].\n"
        "\nParameters\n----------\n"
        "xi : numpy.ndarray, shape (N,)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N,)\n"
        "    Array containing the fluid velocity v\n"
        "idx_wall : integer\n"
        "    Index of the first element outside the bubble\n"
        "idx_shock : integer, optional (default: None)\n"
        "    Index of the first element outside the shock front (use None if"
        " there is no shock front)"
    );

    m.def(
        "compute_energy_density_profile",
        [](py::array_t<double> w, size_t idx_wall, double epsilon) { 
            return vec2np(compute_energy_density_profile(
                np2vec(w), idx_wall, epsilon
            )); 
        },
        py::arg("w"), py::arg("idx_wall"), py::arg("epsilon"),
        "Compute the energy density profile.\n\n"
        "We assume the Bag equation of state, i.e. e = 3/4 w + epsilon.\n"
        "\nParameters\n----------\n"
        "w : numpy.ndarray, shape (N,)\n"
        "    Array containing the enthalpy density w\n"
        "idx_wall : integer\n"
        "    Index of the first element outside the bubble\n"
        "epsilon : float\n    Bag constant\n"
        "\nReturns\n-------\n"
        "e : numpy.ndarray, shape (N,)\n    Energy density profile e(xi)"
    );

    m.def(
        "compute_pressure_density_profile",
        [](py::array_t<double> w, size_t idx_wall, double epsilon) { 
            return vec2np(compute_pressure_density_profile(
                np2vec(w), idx_wall, epsilon
            )); 
        },
        py::arg("w"), py::arg("idx_wall"), py::arg("epsilon"),
        "Compute the pressure density profile.\n\n"
        "We assume the Bag equation of state, i.e. p = 1/4 w - epsilon.\n"
        "\nParameters\n----------\n"
        "w : numpy.ndarray, shape (N,)\n"
        "    Array containing the enthalpy density w\n"
        "idx_wall : integer\n"
        "    Index of the first element outside the bubble\n"
        "epsilon : float\n    Bag constant\n"
        "\nReturns\n-------\n"
        "p : numpy.ndarray, shape (N,)\n    Pressure density profile p(xi)"
    );

    m.def(
        "compute_epsilon",
        &compute_epsilon,
        py::arg("w"), py::arg("alpha"),
        "Compute the Bag constant (vacuum energy density)"
        " epsilon = 3/4 alpha+ w+ = 3/4 alphaN wN.\n"
        "\nParameters\n----------\n"
        "w : float\n    Enthalpy density, either w+ (directly in front) or wN"
        " (far in front of the bubble wall)\n"
        "alpha : float\n    Transition strength, either alpha+ (directly in"
        " front) or alphaN (far in front of the bubble wall)\n"
        "\nReturns\n-------\n"
        "epsilon : float\n    Bag constant epsilon"
    );

    m.def(
        "compute_alpha_N",
        &compute_alpha_N,
        py::arg("sol_type"), 
        py::arg("alpha_plus"), 
        py::arg("xi_wall"), 
        py::arg("step_size"),
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the transition strength alphaN from the value in front of the"
        " bubble wall alpha+.\n"
        "\nParameters\n----------\n"
        "sol_type : HydroSolutionType\n    Type of the solution\n"
        "alpha_plus : float\n    Transition strength alpha+ immediately in"
        " front of the bubble wall\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "step_size : float\n    Step size for the fluid velocity (or"
        " coordinate) when solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "alpha_N : float\n"
        "    Transition strength alphaN far in front of the bubble wall"
    );

    m.def(
        "find_alpha_plus",
        py::overload_cast<
            double, double, double, ProfileSolveMode
        >(&find_alpha_plus),
        py::arg("alpha_N"),
        py::arg("xi_wall"),
        py::arg("step_size"),
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Find the transition strength alpha+ from the strength alphaN in the "
        "false vacuum phase.\n\n"
        "If the solution type is not provided, the function calls"
        " classify_wave_fluid_frame internally.\n"
        "\nParameters\n----------\n"
        "alpha_N : float\n"
        "    Transition strength alphaN far in front of the wall\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "step_size : float\n    Step size for the fluid velocity (or"
        " coordinate) when solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "alpha_plus : float\n    Value of alpha+"
    );

    m.def(
        "find_alpha_plus",
        py::overload_cast<
            HydroSolutionType, 
            double, 
            double, 
            double, 
            ProfileSolveMode
        >(&find_alpha_plus),
        py::arg("sol_type"),
        py::arg("alpha_N"),
        py::arg("xi_wall"),
        py::arg("step_size"),
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Find the transition strength alpha+ from the strength alphaN in the"
        " false vacuum phase.\n"
        "\nParameters\n----------\n"
        "sol_type : HydroSolutionType, optional\n    Type of the solution\n"
        "alpha_N : float\n"
        "    Transition strength alphaN far in front of the wall\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "step_size : float\n    Step size for the fluid velocity (or"
        " coordinate) when solving the differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "alpha_plus : float\n    Value of alpha+"
    );

    m.def(
        "compute_efficiency_factor",
        [](
            py::array_t<double> xi, 
            py::array_t<double> v, 
            py::array_t<double> w, 
            size_t idx_wall, 
            double epsilon
        ) {
            return compute_efficiency_factor(
                np2vec(xi), np2vec(v), np2vec(w), idx_wall, epsilon
            );
        },
        py::arg("xi"), 
        py::arg("v"), 
        py::arg("w"), 
        py::arg("idx_wall"), 
        py::arg("epsilon"),
        "Compute the efficiency factor kappa_sw.\n\n"
        "Compute the efficiency factor kappa_sw for converting vacuum energy to"
        " bulk kinetic energy, Eq. (30) of [Espinosa+ 2010] or Eq. (49) of"
        " [Barni+ 2024].\n"
        "\nParameters\n----------\n"
        "xi : numpy.ndarray, shape (N,)\n"
        "    Array containing the self-similar coordinate xi\n"
        "v : numpy.ndarray, shape (N,)\n"
        "    Array containing the fluid velocity v\n"
        "w : numpy.ndarray, shape (N,)\n"
        "    Array containing the enthalpy density w\n"
        "idx_wall : integer\n"
        "    Index of the first element outside the bubble\n"
        "epsilon : float\n    Bag constant\n"
        "\nReturns\n-------\n"
        "kappa : float\n    Efficiency factor kappa"
    );

    m.def(
        "compute_dofs_ratio_LTE",
        &compute_dofs_ratio_LTE,
        py::arg("alpha_N"), 
        py::arg("xi_wall"), 
        py::arg("step_size")=1e-4,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Compute the ratio Psi of degrees of freedom in local thermal"
        " equilibrium.\n\n"
        "Compute the ration Psi = a-/a+ of degrees of freedom ahead of and"
        " behind the wall. This assumes local thermal equilibrium and uses"
        " entropy conservation to relate the temperatures on both sides of the"
        " wall, i.e. s+ v+ gamma+ = s- v- gamma-.\n"
        "\nParameters\n----------\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "alpha_N : float\n"
        "    Transition strength alphaN (far in front of the wall)\n"
        "step_size : float, optional (default: 1e-4)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "kappa : float\n    Ratio Psi of degrees of freedom across the wall"
    );

    m.def(
        "find_LTE_wall_velocity",
        &find_LTE_wall_velocity,
        py::arg("alpha_N"), 
        py::arg("Psi"), 
        py::arg("xi_min"), 
        py::arg("xi_max"), 
        py::arg("step_size")=1e-4,
        py::arg("solve_mode")=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Find the LTE wall velocity xi_w for a given dof ratio Psi.\n\n"
        "Find the wall velocity xi_w for a given ratio of degrees of freedom"
        " Psi = a-/a+ assuming local thermal equilibrium.\n"
        "\nParameters\n----------\n"
        "alpha_N : float\n"
        "    Transition strength alphaN (far in front of the wall)\n"
        "Psi : float\n    Ratio of degrees of freedom Psi = a-/a+\n"
        "xi_min : float\n    Minimal wall velocity xi_w_min to consider\n"
        "xi_max : float\n    Maximal wall velocity xi_w_max to consider\n"
        "step_size : float, optional (default: 1e-4)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved\n"
        "\nReturns\n-------\n"
        "xi_wall : float\n    Wall velocity xi_w"
    );

    m.def(
        "check_profiles",
        [](
            py::array_t<double> xi, 
            py::array_t<double> profile, 
            bool throw_error
        ) {
            return check_profiles(np2vec(xi), np2vec(profile), throw_error);
        },
        py::arg("xi"), py::arg("profile"), py::arg("throw_error")=false,
        "Check profile for consistency.\n\n"
        "This function returns True if xi and profile have the same size and at"
        " least size 6. Otherwise it throws an error (if throw_error is true)"
        " or returns false.\n"
        "\nParameters\n----------\n"
        "xi :  numpy.ndarray, shape (N,)\n"
        "    Self-similar coordinate xi=r/t of the profile\n"
        "profile : numpy.ndarray, shape (N,)\n"
        "    Vector containing the fluid profile\n"
        "throw_error : bool, optional (default: False)\n"
        "    Whether to throw an error if the test fails\n"
        "\nReturns\n-------\n"
        "valid : bool\nWhether the profile passes the test"
    );
}