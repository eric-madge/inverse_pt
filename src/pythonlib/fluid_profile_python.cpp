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
#include "inverse_pt/fluid_profile.hpp"

namespace py = pybind11;
using namespace inverse_pt;

void bind_fluid_profile(pybind11::module_& m) {

    auto fluid_profile = py::class_<FluidProfile>(
        m, 
        "FluidProfile", 
        "Class representing the fluid profile of a bubble in a cosmological"
        " phase transition."
    );

    fluid_profile.def(
        py::init([](py::array_t<double> xi,
                    py::array_t<double> v,
                    py::array_t<double> w,
                    py::array_t<double> e,
                    py::array_t<double> p,
                    size_t idx_w,
                    py::object idx_sh
        ) {
            return std::make_unique<FluidProfile>( 
                np2vec(xi), 
                np2vec(v), 
                np2vec(w), 
                np2vec(e), 
                np2vec(p), 
                idx_w, 
                obj2idx(idx_sh)
            );
        }),
        py::arg("xi"), 
        py::arg("v"), 
        py::arg("w"), 
        py::arg("e"), 
        py::arg("p"),
        py::arg("idx_w"), 
        py::arg("idx_sh") = 
        py::none(),
        "Constructor for the FluidProfile from precomputed profiles.\n\n"
        "Parameters\n"
        "----------\n"
        "xi : numpy.ndarray, shape (N,)\n"
        "    Self-similar coordinate xi = r/t\n"
        "v : numpy.ndarray, shape (N,)\n"
        "    Fluid velocity profile v(xi)\n"
        "w : numpy.ndarray, shape (N,)\n"
        "    Enthalpy density profile w(xi)\n"
        "e : numpy.ndarray, shape (N,)\n"
        "    Energy density profile e(xi)\n"
        "p : numpy.ndarray, shape (N,)\n"
        "    Pressure density profile p(xi)\n"
        "idx_w : int\n"
        "    Index of the first element outside the bubble\n"
        "idx_sh : int or None, optional\n"
        "    Index of the first element outside the shock front\n"
        "    (use None if there is no shock front)"
    );

    fluid_profile.def(
        py::init<
            double, double, double, profile_calculator::ProfileSolveMode
        >(), 
        py::arg("v_w"), 
        py::arg("alpha_N"), 
        py::arg("step_size") = 1e-3,
        py::arg("solve_mode") 
            = profile_calculator::ProfileSolveMode::D_V_D_XI_FOR_SHOCKS,
        "Constructor for the FluidProfile from the wall velocity xi_w and"
        " strength parameter alpha_N.\n\n"
        "Parameters\n"
        "----------\n"
        "v_w : float\n"
        "    Wall velocity xi_w\n"
        "alpha_N : float\n"
        "    Strength parameter alpha_N (far in front of the wall)\n"
        "step_size : float, optional (default: 0.001)\n"
        "    Step size for the fluid velocity (or coordinate) when solving the"
        " differential equation\n"
        "solve_mode : ProfileSolveMode, optional (default: ProfileSolveMode::"
        "D_V_D_XI_FOR_SHOCKS)\n    ProfileSolveMode switch determining which"
        " differential equation for the velocity profile is solved"
    );

    fluid_profile.def("get_xi_coordinate", 
        [](const FluidProfile& self) {
            return vec2np(self.get_xi_coordinate()); 
        },
        "Get the self-similar coordinate xi=r/t of the profile.\n\n"
        "Returns\n"
        "-------\n"
        "xi : numpy.ndarray\n"
        "    Vector containing the coordinate"
    );

    fluid_profile.def("get_fluid_velocity", 
        [](const FluidProfile& self) {
            return vec2np(self.get_fluid_velocity()); 
        },
        "Get the fluid velocity profile v(xi).\n\n"
        "Returns\n"
        "-------\n"
        "v : numpy.ndarray\n"
        "    Vector containing the velocity"
    );

    fluid_profile.def("get_enthalpy_density", 
        [](const FluidProfile& self) {
            return vec2np(self.get_enthalpy_density()); 
        },
        "Get the enthalpy density profile w(xi).\n\n"
        "Returns\n"
        "-------\n"
        "w : numpy.ndarray\n"
        "    Vector containing the enthalpy"
    );

    fluid_profile.def("get_energy_density", 
        [](const FluidProfile& self) {
            return vec2np(self.get_energy_density()); 
        },
        "Get the energy density profile e(xi).\n\n"
        "Returns\n"
        "-------\n"
        "e : numpy.ndarray\n"
        "    Vector containing the energy"
    );

    fluid_profile.def("get_pressure_density", 
        [](const FluidProfile& self) {
            return vec2np(self.get_pressure_density()); 
        },
        "Get the pressure density profile p(xi).\n\n"
        "Returns\n"
        "-------\n"
        "p : numpy.ndarray\n"
        "    Vector containing the pressure"
    );

    fluid_profile.def("get_wall_index", 
        &FluidProfile::get_wall_index,
        "Get the index of the element in xi corresponding to the bubble wall"
        " position.\n\n"
        "Returns\n"
        "-------\n"
        "idx : int\n"
        "    Index of the bubble wall position (first element outside the wall)"
    );

    fluid_profile.def("get_shock_index",
        [](const FluidProfile& self) -> pybind11::object { 
            return idx2obj(self.get_shock_index()); 
        },
        "Get the index of the element in xi corresponding to the shock front"
        " position.\n\n"
        "Returns\n"
        "-------\n"
        "idx : int\n    Index of the shock front position (first element"
        " outside the shock)"
    );

    fluid_profile.def("get_wall_velocity",
        &FluidProfile::get_wall_velocity,
        "Get the wall velocity xi_w.\n\n"
        "Returns\n"
        "-------\n"
        "xi_w : float\n"
        "    Value of the wall velocity xi_w"
    );

    fluid_profile.def("get_transition_strength",
        &FluidProfile::get_transition_strength,
        "Get the phase transition strength alpha_N.\n\n"
        "Returns\n"
        "-------\n"
        "alpha_N : float\n"
        "    Value of the transition strength alpha_N"
    );

    fluid_profile.def("get_solution_type",
        &FluidProfile::get_solution_type,
        "Get the type of the hydrodynamic solution.\n\n"
        "Returns\n"
        "-------\n"
        "sol_type : HydroSolutionType\n"
        "    Hydrodynamic solution type"
    );

    fluid_profile.def("has_shock",
        &FluidProfile::has_shock,
        "Check whether the profile has a shock front.\n\n"
        "Returns\n"
        "-------\n"
        "b : bool\n"
        "    True if the profile has a shock front, False otherwise"
    );

    fluid_profile.def("get_mean_enthalpy_density",
        &FluidProfile::get_mean_enthalpy_density,
        "Get the mean enthalphy density (far in front of the wall).\n\n"
        "Returns\n"
        "-------\n"
        "w_bar : float\n"
        "    Mean enthalpy density w_bar"
    );

    fluid_profile.def("get_adiabatic_index",
        &FluidProfile::get_adiabatic_index,
        "Get the mean adiabatic index Gamma = w/e in the stable phase.\n\n"
        "Returns\n"
        "-------\n"
        "Gamma : float\n"
        "    Mean adiabatic index Gamma = w/e"
    );

    fluid_profile.def("get_efficiency_factor",
        &FluidProfile::get_efficiency_factor,
        "Get the efficiency factor kappa.\n\n"
        "Returns\n"
        "-------\n"
        "kappa : float\n"
        "    Efficiency factor kappa"
    );

    fluid_profile.def_static("set_verbose",
        &FluidProfile::set_verbose,
        py::arg("b"),
        "Set verbosity.\n\n"
        "Parameters\n"
        "----------\n"
        "b : bool\n"
        "    Boolean value to which the verbosity flag is set"
    );

    fluid_profile.def_static("is_verbose",
        &FluidProfile::is_verbose,
        "Get verbosity.\n\n"
        "Returns\n"
        "-------\n"
        "b : bool\n"
        "    Boolean value to which the verbosity flag is set"
    );
}