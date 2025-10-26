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

#pragma once

#include <pybind11/pybind11.h>

/**
 * @brief Create bindings for settings and constants needed for the profile and
 *     shound shell calculations.
 * @param m Module to which the functions are bound.
 */
void bind_settings(pybind11::module_& m);

/**
 * @brief Create bindings for auxilliary functions needed for the profile and
 * shound shell calculations.
 * @param m Module to which the functions are bound.
 */
void bind_utils(pybind11::module_& m);

/**
 * @brief Create Python pybind11 bindings for the \c FluidProfile class.
 * @param m Module to which the class is bound.
 */
void bind_fluid_profile(pybind11::module_& m);

/**
 * @brief Create Python pybind11 bindings for the \c SoundShellSpectrum class.
 * @param m Module to which the class is bound.
 */
void bind_sound_shell_spectrum(pybind11::module_& m);

/**
 * @brief Create Python pybind11 bindings for the \c profile_calculator
 *     functions.
 * @param m Module to which the functions are bound.
 */
void bind_profile_calculator(pybind11::module_& m);


/**
 * @brief Create Python pybind11 bindings for the \c sound_shell_model
 *     functions.
 * @param m Module to which the functions are bound.
 */
void bind_sound_shell_model(pybind11::module_& m);