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

#include "inverse_pt/utils.hpp"

#include "inverse_pt/settings.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <omp.h>

namespace inverse_pt::utils {

std::vector<double> lin_space(double min_val, double max_val, size_t N) {
    const double step = (max_val - min_val) / (N - 1.0);
    std::vector<double> result(N);
    #pragma omp simd
    for (size_t i=0; i<N; ++i) { result[i] = min_val + i * step; }
    return result;
}

std::vector<double> log_space(double min_val, double max_val, size_t N) {
    const double log_min = std::log(min_val);
    const double log_max = std::log(max_val);
    const double log_step = (log_max - log_min) / (N - 1.0);
    std::vector<double> result(N);
    #pragma omp simd
    for (size_t i=0; i<N; ++i) { result[i] = std::exp(log_min + i * log_step); }
    return result;
}

double integrate_trapez(
    const std::vector<double>& x, const std::vector<double>& y
) {
    size_t N = x.size();
    if ( N != y.size() || N < 1 ) { throw std::invalid_argument(
        "Vectors x and y must have the same size and cannot be empty. Got"
        " sizes " + std::to_string(N) + " and " + std::to_string(y.size()) + "."
    ); }
    double integral = 0.0;
    #pragma omp simd reduction(+:integral)
    for ( size_t i = 0; i < N-1; ++i ) { 
        integral += 0.5 * (y[i+1] + y[i]) * (x[i+1] - x[i]);
    }
    return integral;
}

double integrate_log_trapez(
    const std::vector<double>& x, const std::vector<double>& y
) {
    size_t N = x.size();
    if ( N != y.size() || N < 1 ) { throw std::invalid_argument(
        "Vectors x and y must have the same size and cannot be empty. Got"
        " sizes " + std::to_string(N) + " and " + std::to_string(y.size()) + "."
    ); }
    double integral = 0.0;
    std::vector<double> log_x(N);
    #pragma omp simd
    for ( size_t i = 0; i < N; ++i ) { log_x[i] = std::log(x[i]); }
    #pragma omp simd reduction(+:integral)
    for ( size_t i = 0; i < N-1; ++i ) { 
        integral += 0.5 * (x[i+1]*y[i+1] + x[i]*y[i]) * (log_x[i+1] - log_x[i]);
    }
    return integral;
}

std::pair<double, double> find_peak(
    const std::vector<double>& x, const std::vector<double>& y
) {
    size_t N = x.size();
    if ( N != y.size() || N < 1 ) { throw std::invalid_argument(
        "Vectors x and y must have the same size and cannot be empty. Got" 
        " sizes " + std::to_string(N) + " and " + std::to_string(y.size()) + "."
    ); }

    auto it = std::max_element(y.begin(), y.end());
    size_t index = std::distance(y.begin(), it);
    
    return {x[index], *it};
}

std::vector<double> rescale_vector(
    const std::vector<double>& vector, double factor
) {
    std::vector<double> result(vector.size());
    #pragma omp simd
    for ( size_t i=0; i<vector.size(); ++i ) { result[i] = factor*vector[i]; }
    return result;
}

std::vector<double> rescale_spectrum(
    const std::vector<double>& spectrum, 
    const std::vector<double>& k, 
    int power, 
    double factor
) {
    size_t N = spectrum.size();
    if ( N != k.size() ) {throw std::invalid_argument(
        "The vectors spectrum and k must have the same size. Got " 
        + std::to_string(N) + " and " + std::to_string(k.size()) + "."
    ); }
    std::vector<double> result(N);
    #pragma omp simd
    for ( size_t i=0; i<N; ++i ) {
        result[i] = factor * std::pow(k[i], power) * spectrum[i];
    }
    return result;
}

bool check_spacing(const std::vector<double> x, double tol_abs, double tol_rel)
{
    size_t N = x.size();
    if ( N < 2 ) { return true; }
    double d0 = (x.back()-x.front()) / (N-1.0);
    double tol_fix = std::min(tol_abs, tol_rel*std::abs(d0));
    for (size_t i=0; i<x.size()-1; ++i) {
        double d = x[i+1]-x[i];
        if ( std::abs(d - d0) > std::min(tol_fix, tol_rel*std::abs(d)) ) { 
            return false; 
        }
    }
    return true;
}

void check_failed(const std::string& message) {
    switch(settings::check_mode) {
        case settings::CheckMode::ERROR:
            throw std::runtime_error(message);
        case settings::CheckMode::WARNING:
            std::cerr << "WARNING:\n" << message << std::endl;
            break;
        case settings::CheckMode::SILENT:
        default:
            break;
    }
}


UniformLinearInterpolator::UniformLinearInterpolator(
    const std::vector<double>& x, const std::vector<double>& y
) : 
    x_data(x), 
    y_data(y), 
    size(x.size()), 
    xmin(x.front()), 
    xmax(x.back()), 
    dx((x.back() - x.front()) / (x.size() - 1.0)) 
{
    if (size != y_data.size() || size < 2) {
        throw std::invalid_argument(
            "x and y must have the same size and at least two points."
        );
    }
    if ( ! check_spacing(x) ) {
        throw std::invalid_argument("x must be equally spaced.");
    }
}

double UniformLinearInterpolator::get_xmin() const { return xmin; }
double UniformLinearInterpolator::get_xmax() const { return xmax; }
size_t UniformLinearInterpolator::get_size() const { return size; }
double UniformLinearInterpolator::get_dx() const { return dx; }

UniformCubicHermiteInterpolator::UniformCubicHermiteInterpolator(
    const std::vector<double>& x, 
    const std::vector<double>& y, 
    const std::vector<double>& dy
) : 
    x_data(x), 
    y_data(y), 
    dy_data(dy), 
    size(x.size()), 
    xmin(x.front()), 
    xmax(x.back()), 
    dx((x.back() - x.front()) / (x.size() - 1.0)) 
{
    if (size != y_data.size() || size != dy_data.size() || size < 2) {
        throw std::invalid_argument(
            "x, y and dy must have the same size and at least two points."
        );
    }
    if ( ! check_spacing(x) ) {
        throw std::invalid_argument("x must be equally spaced.");
    }
}

double UniformCubicHermiteInterpolator::get_xmin() const { return xmin; }
double UniformCubicHermiteInterpolator::get_xmax() const { return xmax; }
size_t UniformCubicHermiteInterpolator::get_size() const { return size; }
double UniformCubicHermiteInterpolator::get_dx() const { return dx; }

} // namespace utils