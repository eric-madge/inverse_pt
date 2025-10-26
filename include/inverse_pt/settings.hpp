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

namespace inverse_pt::settings {

    /// report mode options for non-critical consistency checks
    enum class CheckMode{
        ERROR,   ///< throw exception if a check fails
        WARNING, ///< show a warning if a check fails
        SILENT   ///< ignore failed checks
    };

    /// small number corresponding to zero velocity when calculating profiles
    extern double v_zero;
    /// tolerance for considering two values of the coordinate as equal     
    extern double jump_tolerance; 
    /// debug output flag
    extern bool debug;
    /// report mode for non-critical consistency checks
    extern CheckMode check_mode;  

} // namespace settings