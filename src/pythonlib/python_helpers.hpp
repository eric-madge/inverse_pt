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
#include <pybind11/numpy.h>

#include <vector>
#include <stdexcept>

#include "inverse_pt/constants.hpp"

using namespace inverse_pt;

/** 
 * @brief Convert numpy array to std::vector<double>.
 * @param arr 1D numpy array to convert
 * @returns Converted array as std::vector<double>
 */
inline std::vector<double> np2vec(const pybind11::array_t<double>& arr) {
    auto buf = arr.request();
    if (buf.ndim != 1) { throw std::runtime_error("Input array must be 1D."); }
    if (!(arr.flags() & pybind11::array::c_style)) { 
        throw std::runtime_error("Input array must be C-contiguous.");
    }
    const double* ptr = static_cast<const double*>(buf.ptr);
    return std::vector<double>(ptr, ptr + buf.shape[0]);
}

/**
 * @brief Convert std::vector<double> to numpy array.
 * @param vec std::vector<double> to convert
 * @returns Converted vector as 1D numpy array
 */
inline pybind11::array_t<double> vec2np(const std::vector<double>& vec) {
    return pybind11::array_t<double>(vec.size(), vec.data());
}

/**
 * @brief Convert Python object to an index (\c size_t).
 * @details Converts the index (integer) to \c size_t and \c None to
 *     \c constants::NO_INDEX.
 * @param obj Python object (integer) holding the index
 * @returns Index as \c size_t, mapping \c None to \c constants::NO_INDEX
 */
inline size_t obj2idx(pybind11::object obj) {
    return ( 
        obj.is_none() ? constants::NO_INDEX : obj.cast<size_t>() 
    );
}

/**
 * @brief Convert \c size_t to integer index in Python.
 * @details Converse the \c size_t index to an integer index in Python and
 *     \c constants::NO_INDEX to \c None.
 * @param idx \c size_t index in C++
 * @returns Integer index in Python (or \c None if the @p idx is 
 *     \c constants::NO_INDEX)
 */
inline pybind11::object idx2obj(size_t idx) {
    if ( idx == constants::NO_INDEX ) { return pybind11::none(); }
    return pybind11::int_(idx);
}