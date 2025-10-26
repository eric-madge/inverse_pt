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

#include <cstddef>

namespace inverse_pt::constants {
    /// \f$\pi\f$
    inline constexpr double Pi = 3.141592653589793;
    /// speed of sound, \f$c_s = \frac{1}{\sqrt{3}}\f$
    inline constexpr double cs = 0.5773502691896258; 
    /// speed of sound squared, \f$c_s^2 = \frac{1}{3}\f$
    inline constexpr double cs2 = 1.0 / 3.0;         
    /// index value used when not applicable
    /// (e.g. when no shock front is present)
    inline constexpr size_t NO_INDEX = static_cast<size_t>(-1); 
} // namespace constants