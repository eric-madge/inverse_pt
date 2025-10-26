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

#include "inverse_pt_calculator_bindings.hpp"

namespace py = pybind11;

// =============================================
// inverse_pt_calculator module
// =============================================
PYBIND11_MODULE(inverse_pt_calculator, m) {
    m.doc() = "inverse_pt_calculator\n\n"
        "A C++ code for the calculation of the fluid profiles for direct and"
        " inverse cosmological phase transition and the corresponding"
        " graviational wave signal within the sound shell model.\n\n\n"
        "References\n"
        "----------\n"
        "  [Espinosa+ (2010)]\n"
        "    J. R. Espinosa and J. M. No,\n"
        "    \"Energy Budget of Cosmological First-order Phase Transitions,\"\n"
        "    JCAP 06 (2010), 028 [arXiv:1004.4187 [hep-ph]].\n"
        "  [Hindmarsh (2018)]\n" 
        "    M.Hindmarsh,\n"
        "    \"Sound shell model for acoustic gravitational wave production at"
        " a first-order phase transition in the early Universe,\"\n"
        "    Phys. Rev. Lett. 120 (2018) no. 7, 071301 [arXiv: 1608.04735"
        " [astro-ph.CO]].\n"
        "  [Hindmarsh+ (2019)]\n"
        "    M. Hindmarsh and M. Hijazi,\n"
        "    \"Gravitational waves from first order cosmological phase"
        " transitions in the Sound Shell Model,\"\n"
        "    JCAP 12 (2019), 062 [arXiv:1909.10040 [astro-ph.CO]].\n"
        "  [Roper+ (2024)]:\n"
        "    A. Roper Pol, S. Procacci and C. Caprini,\n"
        "    \"Characterization of the gravitational wave spectrum from sound"
        " waves within the sound shell model,\"\n"
        "    Phys. Rev. D 109 (2024) 6, 063531 [arXiv:2308.12943 [gr-qc]].\n"
        "  [Barni+ (2024)]\n"
        "    G. Barni, S. Blasi and M. Vanvlasselaer,\n"
        "    \"The hydrodynamics of inverse phase transitions,\"\n"
        "    JCAP 10 (2024), 042 [arXiv:2406.01596 [hep-ph]].\n"
        "  [Barni+ (2025)]\n"
        "    G. Barni, S. Blasi and M. Vanvlasselaer,\n"
        "    \"Inverse bubbles from broken symmetries,\"\n"
        "    [arXiv:2503.01951 [hep-ph]]."; 

    bind_settings(m);
    bind_utils(m);

    py::module_ pcalc = m.def_submodule(
        "profile_calculator", 
        "Functions to calculate the fluid profile of a bubble in a cosmological"
        " phase transition."
    );
    bind_profile_calculator(pcalc);
    py::module_ ssmodel = m.def_submodule(
        "sound_shell_model", 
        "Functions to calculate the gravitational wave spectrum within the"
        " sound shell model."
    );
    bind_sound_shell_model(ssmodel);

    bind_fluid_profile(m);
    bind_sound_shell_spectrum(m);

} // PYBIND11_MODULE(inverse_pt_calculator, m)
