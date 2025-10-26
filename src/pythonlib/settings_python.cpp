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
#include "inverse_pt/settings.hpp"
#include "inverse_pt/constants.hpp"

namespace py = pybind11;
using namespace inverse_pt;

void bind_settings(pybind11::module_& m) {

    // set/get debug flag for additional command line prints
    m.def("get_debug", []() { return settings::debug; }, "Get the debug flag");
    m.def(
        "set_debug", 
        [](bool b) { settings::debug = b; }, 
        py::arg("b"), 
        "Set the debug flag"
    );

    py::enum_<settings::CheckMode>(
        m, 
        "CheckMode", 
        py::arithmetic(), 
        "Report mode options for non-critical consistency checks"
    )
        .value(
            "ERROR", 
            settings::CheckMode::ERROR, 
            "throw exception if a check fails"
        )
        .value(
            "WARNING", 
            settings::CheckMode::WARNING, 
            "show a warning if a check fails"
        )
        .value(
            "SILENT", 
            settings::CheckMode::SILENT, 
            "ignore failed checks"
        )
    .export_values();
    m.def(
        "get_check_mode", 
        []() { return settings::check_mode; }, 
        "Get the report mode for non-critical consistency checks"
    );
    m.def(
        "set_check_mode", 
        [](settings::CheckMode mode) { settings::check_mode = mode; }, 
        py::arg("mode"), 
        "Set the report mode for non-critical consistency checks"
    );

    // speed of sound
    m.attr("cs")  = py::float_(constants::cs);
    m.attr("cs2") = py::float_(constants::cs2);

}