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
#include "inverse_pt/utils.hpp"

namespace py = pybind11;
using namespace inverse_pt;

void bind_utils(pybind11::module_& m) {

    // linear interpolator class
    auto lin_interp = py::class_<utils::UniformLinearInterpolator>(
        m, 
        "UniformLinearInterpolator",
        "Fast linear interpolator for uniformly spaced data.\n\n"
        "This class performs linear interpolation over a pair of x and y"
        " vectors.\n"
        "**It assumes that the x values are uniformly spaced.**\n"
        "If a query point lies outside the x range, the interpolator returns 0."
    );
    lin_interp.def(
        py::init([](py::array_t<double> x, py::array_t<double> y) {
            return std::make_unique<utils::UniformLinearInterpolator>(
                np2vec(x), np2vec(y)
            );
        }), 
        py::arg("x"), py::arg("y"),
        "Construct a new UniformLinearInterpolator object.\n"
        "\nParameters\n----------\n"
        "x : numpy.ndarray, shape (N,)\n"
        "    Vector of x-values (must be sorted and equally spaced, size >= 2)"
        "\ny : numpy.ndarray, shape (N,)\n"
        "    Vector of y-values corresponding to x-values"
    );
    lin_interp.def(
        "__call__",
        &utils::UniformLinearInterpolator::operator(),
        py::arg("x"),
        "Interpolates the value at a given x position.\n"
        "\nParameters\n----------\n"
        "x : float\n    The x-value to interpolate.\n"
        "\nReturns\n-------\n"
        "y : float\n    Interpolated y-value at x, or 0 if x is out of bounds."
    );
    lin_interp.def_property_readonly(
        "xmin", 
        &utils::UniformLinearInterpolator::get_xmin, 
        "Get minimal x-value."
    );
    lin_interp.def_property_readonly(
        "xmax", 
        &utils::UniformLinearInterpolator::get_xmax, 
        "Get maximal x-value."
    );
    lin_interp.def_property_readonly(
        "size", 
        &utils::UniformLinearInterpolator::get_size, 
        "Get size of the data."
    );
    lin_interp.def_property_readonly(
        "dx", 
        &utils::UniformLinearInterpolator::get_dx, 
        "Get spacing between x-values."
    );

}