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

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <omp.h>

namespace inverse_pt::utils {

/**
 * @brief Linearly spaced vector.
 * 
 * Creates a vector with @p N steps that is linearly spaced between @p min_val
 * and @p max_val.
 * 
 * @param min_val Minimal value
 * @param max_val Maximal value
 * @param N Number of steps
 * @return linearly spaced vector
 */
std::vector<double> lin_space(double min_val, double max_val, size_t N);
    
/**
 * @brief Logarithmically spaced vector.
 * 
 * Creates a vector with @p N steps that is logarithmically spaced between
 * @p min_val and @p max_val.
 * 
 * @param min_val Minimal value
 * @param max_val Maximal value
 * @param N Number of steps
 * @return logarithmically spaced vector
 */
std::vector<double> log_space(double min_val, double max_val, size_t N);

/**
 * @brief Integration using the trapezoidal rule.
 * 
 * This function calculates the integral \f$\int f(x) dx\f$ for the given values
 * of \f$x\f$ and \f$y = f(x)\f$ using the trapezoidal rule
 * \f[
 *  \int\limits_{x_0}^{x_N} dx\,f(x) = \sum\limits_{i=1}^{N-1} \frac{f(x_{i+1})
 *  + f(x_i)}{2} (x_{i_1} - x_i)\,.
 * \f]
 * The function assumes that @p x is ordered.
 * 
 * @param x Grid points \f$x\f$ at which the function is evaluated
 * @param y Function values \f$y = f(x)\f$
 * @result integral \f$\int f(x) dx\f$
 */
double integrate_trapez(
    const std::vector<double>& x, const std::vector<double>& y
);

/**
 * @brief Integration using the trapezoidal rule for logarithmically spaced
 *     \f$x\f$-values.
 * 
 * This function calculates the integral \f$\int f(x) dx\f$ for the given values
 * of \f$x\f$ and \f$y = f(x)\f$ using the trapezoidal rule
 * \f[
 *  \int\limits_{x_0}^{x_N} dx\,f(x)
 *  \int\limits_{\log x_0}^{\log x_N} d(\log x) \,x f(x) 
 *   = \sum\limits_{i=1}^{N-1} \frac{x_{i+1} f(x_{i+1}) +
 *   x_i f(x_i)}{2} (\log x_{i_1} - \log x_i)\,.
 * \f]
 * The function assumes that @p x is ordered and positive.
 * 
 * @param x Grid points \f$x\f$ at which the function is evaluated
 * @param y Function values \f$y = f(x)\f$
 * @result integral \f$\int x f(x) d(\log x)\f$
 */
double integrate_log_trapez(
    const std::vector<double>& x, const std::vector<double>& y
);

/**
 * @brief Find the position and value of the maximum of a discretely sampled
 *     function.
 * 
 * This function searches for the maximum value of a function \f$f(x)\f$, given
 * by the discrete data points \f$(x_i, y_i)\f$, and returns both the
 * corresponding \f$x\f$ position and the maximum value.
 * 
 * @param x Grid points \f$x\f$ at which the function is evaluated
 * @param y Function values \f$y = f(x)\f$ corresponding to the grid points
 * @return A pair \f$(x_{\text{max}}, y_{\text{max}})\f$ where
 *     \f$y_{\text{max}}\f$ is the maximum value of \f$y\f$ and
 *     \f$x_{\text{max}}\f$ is the corresponding grid point
 */
std::pair<double, double> find_peak(
    const std::vector<double>& x, const std::vector<double>& y
);

/** 
 * @brief Rescale (multiply) all entries of a vector by a constant factor.
 * @param vector Vector to rescale
 * @param factor Factor by which the vector is multiplied
 * @return Rescaled vector
 */
std::vector<double> rescale_vector(
    const std::vector<double>& vector, double factor
);

/** 
 * @brief Rescale all entries of the @p spectrum by a constant factor @p factor
 *     and some @p power of the momentum @p k.
 * @param spectrum Spectrum to rescale
 * @param k Momentum
 * @param power Power of the momentum by which the spectrum is multiplied
 * @param factor Constant factor by which the spectrum is multiplied
 *     (optional, default: 1)
 * @return Rescaled vector
 */
std::vector<double> rescale_spectrum(
    const std::vector<double>& spectrum, 
    const std::vector<double>& k, 
    int power, 
    double factor=1.0
);

/**
 * @brief Check wether the vector @p x is evenly spaced.
 * @param x Input vector
 * @param tol_abs Absolute tolerance (optional, default: 1e-12)
 * @param tol_rel Relative tolerance (optional, default: 1e-9)
 */
bool check_spacing(
    const std::vector<double> x, double tol_abs=1e-12, double tol_rel=1e-9
);

/** 
 * @brief Handle failed non-critical consistency checks.
 * 
 * Depending on the value of @ref settings::check_mode, this function throws a
 * \c std::runtime_error (for settings::CheckMode::ERROR), displays a warning
 * (for settings::CheckMode::WARNING), or ignores the failed check
 * (for settings::CheckMode::SILENT).
 * 
 * @param message Error or warning message to display (this has no effect if
 *     @ref settings::check_mode is set to for settings::CheckMode::SILENT)
 */
void check_failed(const std::string& message);

/**
 * @brief Fast linear interpolator for uniformly spaced data.
 * 
 * This class performs linear interpolation over a pair of x and y vectors.
 * **It assumes that the x values are uniformly spaced.**
 * If a query point lies outside the x range, the interpolator returns 0.
 */
class UniformLinearInterpolator {
public:

    /**
     * @brief Construct a new \c UniformLinearInterpolator object.
     * 
     * @param x Vector of x-values (must be sorted and equally spaced,
     *     size >= 2)
     * @param y Vector of y-values corresponding to x-values
     */
    UniformLinearInterpolator(
        const std::vector<double>& x,
        const std::vector<double>& y
    );

    /**
     * @brief Interpolates the value at a given x position.
     * 
     * @param x The x-value to interpolate.
     * @return Interpolated y-value at x, or 0 if x is out of bounds.
     */
    #pragma omp declare simd
    inline double operator()(double x) const {
        const bool out_of_bounds = (x < xmin || x > xmax);

        size_t idx = static_cast<size_t>((x - xmin) / dx);
        idx = std::clamp<size_t>(idx, 0, size - 2);
        if (idx > 0 && x < x_data[idx]) { --idx; }
        else if (idx < size - 2 && x >= x_data[idx + 1]) { ++idx; }

        const double x1 = x_data[idx];
        const double x2 = x_data[idx + 1];
        const double y1 = y_data[idx];
        const double y2 = y_data[idx + 1];

        double result = (y1 * (x2 - x) + y2 * (x - x1)) / (x2 - x1);

        return out_of_bounds ? 0. : result;
    }

    /// Get minimal x-value.
    double get_xmin() const;
    /// Get maximal x-value.
    double get_xmax() const;
    /// Get size of the data.
    size_t get_size() const;
    /// Get spacing between x-values.
    double get_dx() const;

private:
    const std::vector<double> x_data;  ///< x-values of data points
    const std::vector<double> y_data;  ///< y-values of data points
    const size_t size;                 ///< Number of data points
    const double xmin;                 ///< First x-value
    const double xmax;                 ///< Last x-value
    const double dx;                   ///< Average spacing between x-values
};

/**
 * @brief Fast cubic Hermite interpolator for uniformly spaced data.
 * 
 * This class performs cubic Hermite interpolation over a pair of x and y
 * vectors. **It assumes that the x values are uniformly spaced.**
 * If a query point lies outside the x range, the interpolator returns 0.
 */
class UniformCubicHermiteInterpolator {
public:

    /**
     * @brief Construct a new \c UniformCubicHermiteInterpolator object.
     * 
     * @param x Vector of x-values (must be sorted and equally spaced,
     *     size >= 2)
     * @param y Vector of y-values corresponding to x-values
     * @param dy Vector of derivatives at the x-values
     */
    UniformCubicHermiteInterpolator(
        const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& dy
    );

    /**
     * @brief Interpolates the value at a given x position.
     * 
     * @param x The x-value to interpolate.
     * @return Interpolated y-value at x, or 0 if x is out of bounds.
     */
    #pragma omp declare simd
    inline double operator()(double x) const {
        bool out_of_bounds = (x < xmin || x > xmax);

        size_t idx = static_cast<size_t>((x - xmin) / dx);
        idx = std::clamp<size_t>(idx, 0, size - 2);
        if (idx > 0 && x < x_data[idx]) { --idx; }
        else if (idx < size - 2 && x >= x_data[idx + 1]) { ++idx; }

        double t = (x - x_data[idx]) / dx;
        double t2 = t*t;
        double t3 = t*t2;
        double h00 =  2.0*t3 - 3.0*t2 + 1.0;
        double h10 =      t3 - 2.0*t2 + t;
        double h01 = -2.0*t3 + 3.0*t2;
        double h11 =      t3 -     t2;

        return out_of_bounds ? 0. : (
            h00*y_data[idx] + h01 * y_data[idx+1] 
            + dx * ( h10 * dy_data[idx] + h11 * dy_data[idx+1] )
        );
    }

    /// Get minimal x-value.
    double get_xmin() const;
    /// Get maximal x-value.
    double get_xmax() const;
    /// Get size of the data.
    size_t get_size() const;
    /// Get spacing between x-values.
    double get_dx() const;

private:
    const std::vector<double> x_data;  ///< x-values of data points
    const std::vector<double> y_data;  ///< y-values of data points
    const std::vector<double> dy_data; ///< derivative at the data points
    const size_t size;                 ///< Number of data points
    const double xmin;                 ///< First x-value
    const double xmax;                 ///< Last x-value
    const double dx;                   ///< Average spacing between x-values
};

} // namespace utils