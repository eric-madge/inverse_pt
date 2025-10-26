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

#include "inverse_pt/profile_calculator.hpp"

#include "inverse_pt/constants.hpp"
#include "inverse_pt/settings.hpp"
#include "inverse_pt/utils.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/math/tools/roots.hpp>
#include <omp.h>

using namespace inverse_pt::constants;
using namespace inverse_pt::settings;

/// relative tolerance when finding alpha_+ from alpha_N
constexpr double alpha_rel_tol = 1e-6; 
/// absolute tolerance when finding alpha_+ from alpha_N
constexpr double alpha_abs_tol = 1e-9; 
/// relative tolerance when finding wall velocity from LTE condition
constexpr double xi_wall_rel_tol = 1e-6;
/// absolute tolerance when finding wall velocity from LTE condition 
constexpr double xi_wall_abs_tol = 1e-9;

namespace inverse_pt::profile_calculator {

const std::string to_string(HydroSolutionType type) {
    switch (type) {
        case HydroSolutionType::Deflagration:        
            return "Deflagration";
        case HydroSolutionType::Detonation:          
            return "Detonation";
        case HydroSolutionType::Hybrid:              
            return "Hybrid";
        case HydroSolutionType::InverseDeflagration: 
            return "Inverse Deflagration";
        case HydroSolutionType::InverseDetonation:   
            return "Inverse Detonation";
        case HydroSolutionType::InverseHybrid:       
            return "Inverse Hybrid";
        default:                                     
            return "Unknown";
    }
}

std::ostream& operator<<(std::ostream& os, HydroSolutionType type) { 
    return os << to_string(type); 
}

bool are_equal(double x, double y, double eps_rel, double eps_abs) {
    double precision = std::max(
        eps_abs, 
        eps_rel * std::max(std::abs(x), std::abs(y))
    );
    return std::abs(x - y) < precision;
}

double gamma(double v) { return 1.0 / std::sqrt(1.0 - v * v); }
double mu(double xi, double v) { return (xi - v) / (1.0 - v * xi); }
bool check_velocity(double v) { return (v > 0.0 && v < 1.0); }
double shock_condition(double xi, double v) { return mu(xi, v) * xi - cs2; }
double jouguet_condition(double xi, double v) {return mu(xi, v)*mu(xi, v)-cs2;}

std::pair<double, double> get_v_plus_from_v_minus(double vm, double alpha) {
    double term = 0.5 * vm + 1.0 / (6.0 * vm);
    double disc = std::sqrt(
        term * term + alpha * alpha + 2.0 / 3.0 * alpha - 1.0 / 3.0
    );
    double denom = 1.0 + alpha;
    double vplus1 = (term - disc) / denom;
    double vplus2 = (term + disc) / denom;
    // handle limits alpha -> -1: '-' branch is finite, '+' branch is infinite
    if ( ! std::isfinite(1.0 / denom) ) { 
        return {2.0 / 3.0 / term, 2.0 * term / denom - 2.0 / 3.0 / term}; 
    }
    return {vplus1, vplus2};
}
std::pair<double, double> get_v_minus_from_v_plus(double vp, double alpha) {
    double term = 0.5 * (1. + alpha) * vp + (1.0 - 3.0 * alpha) / (6. * vp);
    double disc = std::sqrt(term * term - 1.0 / 3.0);
    double vminus1 = term - disc;
    double vminus2 = term + disc;
    return {vminus1, vminus2};
}

double compute_alpha_plus(double vp, double vm) {
    double factor1 = (3.0 * vp * vm - 1.0);
    double factor2 = (vp - vm);
    double denominator = 3.0 * (1.0 - vp * vp) * vm;
    return factor1 * factor2 / denominator;
}

double direct_jouguet_velocity(double alpha) {
    if ( alpha < 0.0 ) { throw std::invalid_argument(
        "Attempted to calculate Jouguet velocity for negative alpha: " + 
        std::to_string(alpha) + "."
    ); }
    auto [v1, v2] = get_v_plus_from_v_minus(cs, alpha);
    // v_J > c_s -> use '+' branch
    if ( v2 < cs || !check_velocity(v2) ) { throw std::runtime_error(
        "Invalid value for alpha in calculation of direct Jouguet velocity: " + 
        std::to_string(alpha) + 
        ". Got v_J = " + std::to_string(v2) + "."
    ); }
    return v2; 
}
double inverse_jouguet_velocity(double alpha) {
    if ( alpha > 0.0 ) { throw std::invalid_argument(
        "Attempted to calculate inverse Jouguet velocity for positive alpha: " + 
        std::to_string(alpha) + "."
    ); }
    auto [v1, v2] = get_v_minus_from_v_plus(cs, alpha);
    // inverse v_J < c_s -> use '-' branch
    if ( v1 > cs || ! check_velocity(v1) ) { throw std::runtime_error(
        "Invalid value for alpha in calculation of inverse Jouguet velocity: " + 
        std::to_string(alpha) + 
        ". Got v_J = " + std::to_string(v1) + "."
    ); }
    return v1; 
}

double alpha_max_detonation(
    double xi_w, 
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    if ( xi_w < cs ) { throw std::invalid_argument(
        "The wall velocity for direct detonations must be supersonic (xi_w >"
        " cs). Got " + std::to_string(xi_w) + "."
    ); }
    return compute_alpha_plus(xi_w, cs);
}

double alpha_max_deflagration(
    double xi_w, 
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    if ( xi_w > cs ) { throw std::invalid_argument(
        "The wall velocity for direct deflagrations must be subsonic (xi_w <"
        " cs). Got " + std::to_string(xi_w) + "."
    ); }
    if ( return_alpha_N ) {
        return compute_alpha_N(
            HydroSolutionType::Deflagration, 
            cs2-alpha_abs_tol, 
            xi_w, 
            step_size,
            solve_mode
        );
    }
    return cs2;
}

double alpha_max_hybrid(
    double xi_w, 
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    if ( xi_w < cs ) { throw std::invalid_argument(
        "The wall velocity for direct hybrids must be supersonic (xi_w > cs)."
        " Got " + std::to_string(xi_w) + "."
    ); }
    if ( return_alpha_N ) {
        return compute_alpha_N(
            HydroSolutionType::Hybrid, 
            cs2-alpha_abs_tol,
            xi_w, 
            step_size,
            solve_mode
        );
    }
    return cs2;
}

double alpha_min_hybrid(
    double xi_w, 
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    if ( xi_w < cs ) { throw std::invalid_argument(
        "The wall velocity for direct hybrids must be supersonic (xi_w > cs)."
        " Got " + std::to_string(xi_w) + "."
    ); }
    double alpha_min = cs2 * (xi_w - cs)*(xi_w - cs) / (xi_w*xi_w - cs2*cs2);
    if ( return_alpha_N ) {
        return compute_alpha_N(
            HydroSolutionType::Hybrid, 
            alpha_min+alpha_abs_tol, 
            xi_w, 
            step_size,
            solve_mode
        );
    }
    return alpha_min;
}

double alpha_min_inverse_deflagration(
    double xi_w, 
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    if ( xi_w < cs ) { throw std::invalid_argument(
        "The wall velocity for inverse deflagrations must be supersonic (xi_w >"
        " cs). Got " + std::to_string(xi_w) + "."
    ); }
    return  (cs2 - xi_w) / (1.0 + xi_w);
}

double alpha_min_inverse_hybrid(
    double xi_w, 
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    if ( xi_w > cs ) { throw std::invalid_argument(
        "The wall velocity for inverse hybrids must be subsonic (xi_w < cs)."
        " Got " + std::to_string(xi_w) + "."
    ); }
    if ( xi_w < cs2 ) { throw std::invalid_argument(
        "The wall velocity for inverse hybrids must larger than \f$c_s^2\f$."
        " Got " + std::to_string(xi_w) + "."
    ); }
    double alpha_min = 1.0 - 2.0 * cs;
    if ( return_alpha_N ) {
        return compute_alpha_N(
            HydroSolutionType::InverseHybrid, 
            alpha_min+alpha_abs_tol, 
            xi_w, 
            step_size,
            solve_mode
        );
    }
    return alpha_min;
}

double alpha_min_inverse_detonation(
    double xi_w, 
    bool return_alpha_N,
    double step_size,
    ProfileSolveMode solve_mode
) {
    if ( xi_w > cs ) { throw std::invalid_argument(
        "The wall velocity for inverse hybrids must be subsonic (xi_w < cs)."
        " Got " + std::to_string(xi_w) + "."
    ); }
    double alpha_min = compute_alpha_plus(cs, xi_w);
    if ( return_alpha_N ) {
        return compute_alpha_N(
            HydroSolutionType::InverseDetonation, 
            alpha_min+alpha_abs_tol, 
            xi_w, 
            step_size,
            solve_mode
        );
    }
    return alpha_min;
}

double alpha_max(
    HydroSolutionType sol_type, 
    double xi_w, 
    bool return_alpha_N, 
    double step_size,
    ProfileSolveMode solve_mode
) {
    switch ( sol_type ) {
        case HydroSolutionType::Detonation:
            return alpha_max_detonation(
                xi_w, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::Deflagration: 
            return alpha_max_deflagration(
                xi_w, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::Hybrid: 
            return alpha_max_hybrid(
                xi_w, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::InverseDeflagration:
        case HydroSolutionType::InverseDetonation:
            return 0.0;
        case HydroSolutionType::InverseHybrid: 
            return alpha_min_inverse_detonation(
                xi_w, return_alpha_N, step_size, solve_mode
            );
        default:
            throw std::runtime_error(
                "Unknown solution type: " + to_string(sol_type) +"."
            );
    }
}

double alpha_min(
    HydroSolutionType sol_type, 
    double xi_w, 
    bool return_alpha_N, 
    double step_size,
    ProfileSolveMode solve_mode
) {
    switch ( sol_type ) {
        case HydroSolutionType::Detonation: 
        case HydroSolutionType::Deflagration: 
            return 0.0;
        case HydroSolutionType::Hybrid: 
            return alpha_min_hybrid(
                xi_w, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::InverseDeflagration:
            return alpha_min_inverse_deflagration(
                xi_w, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::InverseDetonation:
            return alpha_min_inverse_detonation(
                xi_w, return_alpha_N, step_size, solve_mode
            );
        case HydroSolutionType::InverseHybrid: 
            return alpha_min_inverse_hybrid(
                xi_w, return_alpha_N, step_size, solve_mode
            );
        default:
            throw std::runtime_error(
                "Unknown solution type: " + to_string(sol_type) +"."
            );
    }
}

double alpha_max(
    double xi_w, 
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    if ( xi_w < cs ) { 
        return alpha_max_deflagration(
            xi_w, return_alpha_N, step_size, solve_mode
        ); 
    }
    if ( return_alpha_N || xi_w < direct_jouguet_velocity(cs2) ) { 
        return alpha_max_hybrid(xi_w, return_alpha_N, step_size, solve_mode); 
    }
    return alpha_max_detonation(xi_w, return_alpha_N, step_size, solve_mode);
}

double alpha_min(
    double xi_w, 
    bool return_alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    if ( xi_w > cs ) { 
        return alpha_min_inverse_deflagration(
            xi_w, return_alpha_N, step_size, solve_mode
        ); 
    }
    if ( xi_w > cs2 ) { 
        return alpha_min_inverse_hybrid(
            xi_w, return_alpha_N, step_size, solve_mode
        ); 
    }
    return alpha_min_inverse_detonation(
        xi_w, return_alpha_N, step_size, solve_mode
    );
}

HydroSolutionType classify_wave_fluid_frame(
    double xi_wall, 
    double alpha, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    // direct solutions 
    if ( alpha >= 0.0 ) {
        // detonations can be identified directly since alphaN = alpha+
        if ( xi_wall > cs && alpha < alpha_max_detonation(
            xi_wall, true, step_size, solve_mode
        ) ) { return HydroSolutionType::Detonation; }
        // for deflagrations: check if the value of alphaN is consistent
        if ( xi_wall < cs && alpha < alpha_max_deflagration(
            xi_wall, true, step_size, solve_mode
        ) ) { return HydroSolutionType::Deflagration; }
        // for hybrids: check if the value of alphaN is consistent
        // no need to check if consistent with alpha_min_hybrid, as this should
        // coincide with alpha_max_deflagration for alpha_N
        if ( xi_wall > cs && alpha < alpha_max_hybrid(
            xi_wall, true, step_size, solve_mode
        ) ) { return HydroSolutionType::Hybrid; }
    } 
    // inverse solutions (alpha < 0)
    // inverse deflagrations can be identified directly since alphaN = alpha+
    if ( xi_wall > cs && alpha > alpha_min_inverse_deflagration(
        xi_wall,true, step_size, solve_mode
    ) ){ return HydroSolutionType::InverseDeflagration; } 
    // for inverse detonations: check that alphaN is consistent
    // (this then automatically satisfies xi_wall < xiJ_inv)
    if ( xi_wall < cs && alpha > alpha_min_inverse_detonation(
        xi_wall, true, step_size, solve_mode
    ) ) { return HydroSolutionType::InverseDetonation; }
    // for inverse hybrids: check if the value of alphaN is consistent
    if ( xi_wall < cs && alpha > alpha_min_inverse_hybrid(
        xi_wall, true, step_size, solve_mode
    ) ) { return HydroSolutionType::InverseHybrid; }
    throw std::invalid_argument(
        "The parameters xi_wall = " + std::to_string(xi_wall) + 
        " and alpha_N = " + std::to_string(alpha) +
        " do not correspond to a valid hydrodynamic solution."
    );
}

HydroSolutionType classify_wave_wall_frame(double vp, double vm) {
    if ( check_velocity(vp) && check_velocity(vm) ) {
        if (vp < vm) {
            if ( vm < cs ) { return HydroSolutionType::Deflagration; }
            if ( vp > cs ) { return HydroSolutionType::InverseDeflagration; }
            // These are actually strong deflagrations, but can be turned into
            // hybrid by setting v_- to c_s:
            if ( vp * vm < cs2 ) { return HydroSolutionType::Hybrid; } 
            // These are actually strong inverse deflagrations, but can be
            // turned into inverse hybrid by setting v_+ to c_s:
            if ( vp * vm > cs2 ) { return HydroSolutionType::InverseHybrid; }
        } 
        // if ( vp >= vm )
        if ( vm > cs ) { return HydroSolutionType::Detonation; }
        if ( vp < cs ) { return HydroSolutionType::InverseDetonation; }
    }
    throw std::runtime_error(
        "Could not identify solution type. Invalid velocities: " 
        + std::to_string(vp) + " and " + std::to_string(vm) + "."
    );
}

std::pair<double, double> compute_fluid_velocities(
    double xi_wall, double alpha, HydroSolutionType sol_type
) {
    double v_plus, v_minus;
    std::pair<double, double> v_branches;
    if ( debug ) { 
        std::cout << "[DEBUG: compute_fluid_velocities] (xi_w, alpha_+) = ("
                  << xi_wall << ", " << alpha << ") -> "; 
    }
    switch( sol_type ) {
        // deflagration: v_ = xi_w < cs, v_+ < cs ('-' branch)
        case HydroSolutionType::Deflagration: 
            v_minus = xi_wall;
            v_branches = get_v_plus_from_v_minus(v_minus, alpha);
            v_plus = v_branches.first;
            if ( debug ) { std::cout << "Deflagration (xi_w < c_s)"; }
            break;
        // hybrid: v_- = c_s, v_+ < cs ('-' branch)
        // note: v_- fixed to cs, return xi_w instead of v_-
        case HydroSolutionType::Hybrid:
            v_minus = xi_wall;
            v_branches = get_v_plus_from_v_minus(cs, alpha);
            v_plus = v_branches.first;
            if ( debug ) { 
                std::cout << "Hybrid (c_s < xi_w < v_J = " 
                          << direct_jouguet_velocity(alpha) << ")"; 
            }
            break;
        // detonation: v_+ = xi_w > xi_J > cs, v_- > cs ('+' branch)
        case HydroSolutionType::Detonation:
            v_plus = xi_wall; 
            v_branches = get_v_minus_from_v_plus(v_plus, alpha);
            v_minus = v_branches.second;
            if ( debug ) {
                std::cout << "Detonation (xi_w > v_J = " 
                          << direct_jouguet_velocity(alpha) << ")";
            }
            break;
        // inverse detonation: v_- = xi_w < v_Jinv < cs, v_+ < c_s ('-' branch)
        case HydroSolutionType::InverseDetonation:
            v_minus = xi_wall;
            v_branches = get_v_plus_from_v_minus(v_minus, alpha);
            v_plus = v_branches.first;
            if ( debug ) { 
                std::cout << "Inverse Detonation (xi_w < v_Jinv = "
                          << inverse_jouguet_velocity(alpha) << ")"; 
            }
            break;
        // inverse hybrid: v_+ = cs, v_- > cs ('+' branch)
        // note: v_+ fixed to cs, return xi_w instead of v_+
        case HydroSolutionType::InverseHybrid:
            v_plus = xi_wall; 
            v_branches = get_v_minus_from_v_plus(cs, alpha);
            v_minus = v_branches.second;
            if ( debug ) { 
                std::cout << "Inverse Hybrid (c_s > xi_w > v_Jinv = " 
                          << inverse_jouguet_velocity(alpha) << ")"; 
            }
            break; 
        // inverse deflagrationL v_+ = xi_w > cs, v_- > cs ('+' branch)
        case HydroSolutionType::InverseDeflagration:
            v_plus = xi_wall; 
            v_branches = get_v_minus_from_v_plus(v_plus, alpha);
            v_minus = v_branches.second;
            if ( debug ) { std::cout << "Inverse Deflagration (xi_w > c_s)"; }
            break; 
        default:
            throw std::runtime_error(
                "Unknown solution type: " + to_string(sol_type) +"."
            );
    }
    if ( debug ) { 
        std::cout << " -> v_+ = " << v_plus << ", v_- = " << v_minus 
                  << std::endl; 
    }
    // check if the velocities are valid 
    // (this should be fine for alpha_+ in the allowed range)
    if ( !(check_velocity(v_plus) && check_velocity(v_minus)) ) {
        throw std::runtime_error(
            "Invalid kinematic configuration: xi_w = " + std::to_string(xi_wall)
            + ", alpha_+ = " + std::to_string(alpha) + ". Found v_+ = "
            + std::to_string(v_plus) + " and v_- = " + std::to_string(v_minus)
            + "."
        );
    }
    // check if classification in terms of v_+ and v_- gives the correct result
    // note: for (inverse) hybrids, classification with xi_w instead of cs for 
    //       v_- (v_+) moves the solution into strong (inverse) deflagration 
    //       regime, but the code identifies these as hybrids
    if ( sol_type != classify_wave_wall_frame(v_plus, v_minus) ) {
        utils::check_failed(
            "[CHECKS: compute_fluid_velocities] "
            "Solution type changed after calculating the fluid velocities for"
            " xi_w = " + std::to_string(xi_wall) + " and alpha_+ = " 
            + std::to_string(alpha)+ ": " + to_string(sol_type) + " vs. " 
            + to_string(classify_wave_wall_frame(v_plus, v_minus)) + "."
        );
    }
    if ( ( 
        sol_type == HydroSolutionType::Hybrid &&
        sol_type != classify_wave_wall_frame(v_plus, cs) 
    ) || ( 
        sol_type == HydroSolutionType::InverseHybrid &&
        sol_type != classify_wave_wall_frame(cs, v_minus)  
    ) ) {
        utils::check_failed(
            "[CHECKS: compute_fluid_velocities] "
            "Hybrid or Inverse Hybrid no longer recongnized if fluid velocity"
            " is set to speed of sound for xi_w = " + std::to_string(xi_wall) + 
            " and alpha_+ = " + std::to_string(alpha)+ ": " 
            + to_string(sol_type) + " vs. " + to_string(
                sol_type == HydroSolutionType::Hybrid ?
                classify_wave_wall_frame(v_plus, cs) : 
                classify_wave_wall_frame(cs, v_minus)
            ) + "."
        );
    }
    return {v_plus, v_minus};
}

double dxi_dv(double v, double xi) {
    double term1 = (xi - v);
    double term2 = (1.0 - xi * v);
    double term3 = (1.0 - v * v);
    double dxi_dt = xi * (term1 * term1 - cs2 * term2 * term2);
    double dv_dt = 2.0 * v * cs2 * term3 * term2;
    return dxi_dt / dv_dt;
}

double dv_dxi(double xi, double v) {
    double term1 = (xi - v);
    double term2 = (1.0 - xi * v);
    double term3 = (1.0 - v * v);
    double dxi_dt = xi * (term1 * term1 - cs2 * term2 * term2);
    double dv_dt = 2.0 * v * cs2 * term3 * term2;
    return  dv_dt / dxi_dt;
}

int integrate_xi_of_v(
    std::vector<double>& xi_vec, 
    std::vector<double>& v_vec, 
    double v_stop, 
    double dv, 
    bool stop_at_shock
) {
    if ( xi_vec.size() < 1 || xi_vec.size() != v_vec.size() ) {
        throw std::runtime_error(
            "The solution vectors must have equal size and contain the initial "
            "condition. Got sizes " + std::to_string(xi_vec.size()) + " and " 
            + std::to_string(v_vec.size()) + "."
        );
    }
    // use the last value stored in xi and v as initial conditions
    // (and apped the results to the vectors)
    double xi = xi_vec.back(), v = v_vec.back();
    // decide whether to evolve in increasing or decreasing velocity direction
    if ( v_stop < v ) { dv = -std::abs(dv); } else { dv = std::abs(dv); }
    // determine when to stop (depending on the integration direction)
    std::function<bool(double)> advance_condition;
    if ( dv > 0.0 ) { 
        advance_condition = [&](double v) { return v <= v_stop; };
    } else {
        advance_condition = [&](double v) { return v >= v_stop; };
    }
    if ( debug ) {
        std::cout << "[DEBUG: integrate_xi_of_v] Integrating " 
                  << (dv > 0.0 ? "forward" : "backward") << " from (v, xi) = (" 
                  << v << ", " << xi << ") to v = " << v_stop << "." 
                  << std::endl;
    }
    double dxi, k1,k2,k3,k4;
    double shock1, shock2;
    // evolve differential equation with 4th-order Runge-Kutta
    while ( advance_condition(v) ) {
        // determine step size in xi
        k1 = dv * dxi_dv(v, xi);
        k2 = dv * dxi_dv(v + 0.5*dv, xi + 0.5*k1);
        k3 = dv * dxi_dv(v + 0.5*dv, xi + 0.5*k2);
        k4 = dv * dxi_dv(v + dv, xi + k3);
        dxi = (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
        // stop if the self-similar coordinate becomes negative or super-luminal
        if ( !check_velocity(xi+dxi) ) { 
            if ( debug ) { 
                std::cout << "[DEBUG: integrate_xi_of_v] Stopped at (v, xi) = "
                          << v << ", " << xi << ") " << "as xi = " << xi+dxi 
                          << " < 0 in next step." << std::endl;
            }
            return -1; 
        }
        // stop if the fluid velocity becomes super-luminal
        if ( std::abs(v+dv) >= 1.0 ) { 
            if ( debug ) { 
                std::cout << "[DEBUG: integrate_xi_of_v] Stopped at (v, xi) = " 
                          << v << ", " << xi << ") " << "as v = " << v+dv 
                          << " becomes super-luminal in next step." 
                          << std::endl;
            }
            return -2; 
        }
        // stop if a shock front is crossed
        if ( stop_at_shock ) {
            shock1 = shock_condition(xi,v);
            shock2 = shock_condition(xi+dxi, v+dv);
            if ( shock1 * shock2 <= 0.0 ) {
                // crossed a shock: determine shock position by linear 
                // interpolation and stop, make sure that we don't overshoot
                // v_stop
                double t = (
                    advance_condition(v + dv) ?
                    shock1 / (shock1 - shock2) :
                    (v_stop - v)/dv
                );
                xi = xi + t * dxi;
                v = v + t * dv;
                xi_vec.emplace_back(xi);
                v_vec.emplace_back(v);
                if ( debug ) {
                    std::cout << "[DEBUG: integrate_xi_of_v] Encountered shock"
                              << " at (v, xi) = (" << v << ", " << xi << ")."
                              << std::endl;
                }
                return 1;
            }
        }
        // step forward
        xi += dxi; v += dv;
        if ( !advance_condition(v) ) { break; }
        xi_vec.emplace_back(xi);
        v_vec.emplace_back(v);
    }

    return 0;
}

int integrate_v_of_xi(
    std::vector<double>& xi_vec, 
    std::vector<double>& v_vec, 
    double v_stop, 
    double dxi, 
    bool stop_at_shock
) {
    if ( xi_vec.size() < 1 || xi_vec.size() != v_vec.size() ) {
        throw std::runtime_error(
            "The solution vectors must have equal size and contain the initial"
            " condition. Got sizes " + std::to_string(xi_vec.size()) + " and " 
            + std::to_string(v_vec.size()) + "."
        );
    }
    // use the last value stored in xi and v as initial conditions
    // (and apped the results to the vectors)
    double xi = xi_vec.back(), v = v_vec.back();
    // decide whether to evolve in increasing or decreasing coordinate direction
    dxi = (v - v_stop) / dv_dxi(xi,v) > 0.0 ? -std::abs(dxi) : std::abs(dxi);
    // determine when to stop (depending on the integration direction)
    std::function<bool(double)> advance_condition;
    if ( v < v_stop ) { 
        advance_condition = [&](double v) { return v <= v_stop; };
    } else {
        advance_condition = [&](double v) { return v >= v_stop; };
    }
    if ( debug ) {
        std::cout << "[DEBUG: integrate_v_of_xi] Integrating " 
                  << (dxi > 0.0 ? "forward" : "backward") << " from (xi, v) = ("
                  << xi << ", " << v << ") to v = " << v_stop << "." 
                  << std::endl;
    }
    double dv, k1,k2,k3,k4;
    double shock1, shock2;
    double jouguet1 = jouguet_condition(xi,v);
    double jouguet2;
    // evolve differential equation with 4th-order Runge-Kutta
    while ( advance_condition(v) ) {
        // determine step size in xi
        k1 = dxi * dv_dxi(xi, v);
        k2 = dxi * dv_dxi(xi + 0.5*dxi, v + 0.5*k1);
        k3 = dxi * dv_dxi(xi + 0.5*dxi, v + 0.5*k2);
        k4 = dxi * dv_dxi(xi + dxi, v + k3);
        dv = (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
        // stop if the self-similar coordinate becomes negative or super-luminal
        if ( !check_velocity(xi+dxi) ) { 
            if ( debug ) { 
                std::cout << "[DEBUG: integrate_v_of_xi] Stopped at (xi, v) = "
                          << xi << ", " << v << ") " << "as xi = " << xi+dxi 
                          << " becomes unphysical in next step." << std::endl;
            }
            return -1; 
        }
        // stop if the fluid velocity becomes super-luminal
        if ( std::abs(v+dv) >= 1.0 ) { 
            if ( debug ) { 
                std::cout << "[DEBUG: integrate_v_of_xi] Stopped at (xi, v) = "
                          << xi << ", " << v << ") " << "as v = " << v+dv 
                          << " becomes super-luminal in next step." 
                          << std::endl;
            }
            return -2; 
        }
        // stop if a shock front is crossed
        if ( stop_at_shock ) {
            shock1 = shock_condition(xi,v);
            shock2 = shock_condition(xi+dxi, v+dv);
            if ( shock1 * shock2 <= 0.0 ) {
                // crossed a shock: determine shock position by linear 
                // interpolation and stop, make sure that we don't overshoot
                // v_stop
                double t = (
                    advance_condition(v + dv) ?
                    shock1 / (shock1 - shock2) :
                    (v_stop - v)/dv
                );
                xi = xi + t * dxi;
                v = v + t * dv;
                xi_vec.emplace_back(xi);
                v_vec.emplace_back(v);
                if ( debug ) {
                    std::cout << "[DEBUG: integrate_v_of_xi] Encountered shock"
                              << " at (xi, vi) = (" << xi << ", " << v << ")." 
                              << std::endl;
                }
                return 1;
            }
        }
        // stop if the Jouguet front is crossed
        // this probably never occurs, the profile will probably evolve along 
        // the front, but anyway we are not going to solve the equations in this
        // direction.
        // jouguet1 = jouguet_condition(xi, v);
        jouguet2 = jouguet_condition(xi+dxi, v+dv);
        if ( jouguet1 * jouguet2 <= 0.0 ) {
            if ( debug ) {
                std::cout << "[DEBUG: integrate_v_of_xi] Stopped at (xi, v) = (" 
                          << xi << ", " << v << ")" << " as the Jouguet front"
                          << " is crossed in the next step." << std::endl;
            }
            return -3;
        }
        // step forward
        xi += dxi; v += dv;
        if ( !advance_condition(v) ) { break; }
        xi_vec.emplace_back(xi);
        v_vec.emplace_back(v);
    }
    return 0;
}

int integrate_velocity_profile(
    std::vector<double>& xi_vec,
    std::vector<double>& v_vec,
    double v_stop,
    double step_size,
    bool stop_at_shock,
    ProfileSolveMode solve_mode
) {
    switch ( solve_mode ) {
        case ProfileSolveMode::D_XI_D_V:
            return integrate_xi_of_v(
                xi_vec, v_vec, v_stop, step_size, stop_at_shock
            );
        case ProfileSolveMode::D_V_D_XI:
            return integrate_v_of_xi(
                xi_vec, v_vec, v_stop, step_size, stop_at_shock
            );
        case ProfileSolveMode::D_V_D_XI_FOR_SHOCKS:
            return (
                stop_at_shock ?
                integrate_v_of_xi(
                    xi_vec, v_vec, v_stop, step_size, stop_at_shock
                ) :
                integrate_xi_of_v(
                    xi_vec, v_vec, v_stop, step_size, stop_at_shock
                )
            );
    }
    throw std::runtime_error(
        "profile_calculator::integrate_velocity_profile got unknown solve mode."
    );
}

size_t compute_velocity_profile_direct_detonation(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double v_minus, 
    double step_size,
    ProfileSolveMode solve_mode
) {
    // check solution type
    if ( 
        classify_wave_wall_frame(v_plus, v_minus) 
        != HydroSolutionType::Detonation 
    ) {
        throw std::runtime_error(
            "The solution with v_+ = " + std::to_string(v_plus) + " and v_- = "
            + std::to_string(v_minus) + " is not a detonation but of type " 
            + to_string(classify_wave_wall_frame(v_plus, v_minus)) + "."
        );
    }
    // calculate fluid velocity in front / behind bubble in the fluid frame
    double xi_wall = v_plus;
    double v0 = mu(xi_wall, v_minus);
    if ( debug ) { 
        std::cout << "[DEBUG: compute_velocity_profile] Direct Detonation:"
                  << " xi_w = " << xi_wall << ", v(xi_w) = " << v0 << "." 
                  << std::endl; 
    }
    
    // solve fluid equations
    int status; // integration return status
    xi.clear(); v.clear();
    // add initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0);
    // integrate from wall to sound speed
    status = integrate_velocity_profile(
        xi, v, v_zero, step_size, false, solve_mode
    );
    if ( status != 0 ) {
        utils::check_failed(
            "[CHECKS: compute_velocity_profile_direct_detonation] "
            "Expected integration return status 0 (reached v_stop), but got "
            + std::to_string(status) + " (" 
            + (status==1 ? "encountered shock" : "entered unphysical region") 
            + ") for v_+ = " + std::to_string(v_plus) + ", v_- = " 
            + std::to_string(v_minus) + "."
        );
    }
    // add first two elemnts
    xi.push_back(cs); xi.push_back(0.);
    v.push_back(0.); v.push_back(0.);
    // revert order 
    std::reverse(xi.begin(), xi.end());
    std::reverse(v.begin(), v.end());
    // set wall velocity index
    size_t idx_wall = xi.size(); // index of first element outside the bubble
    // add last two elements
    xi.push_back(xi.back()); xi.push_back(1.);
    v.push_back(0.); v.push_back(0.);
    
    // return index of wall velocity
    return idx_wall;
}

size_t compute_velocity_profile_direct_deflagration(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double v_minus, 
    double step_size,
    ProfileSolveMode solve_mode
) {
    // check solution type
    if ( 
        classify_wave_wall_frame(v_plus, v_minus) 
        != HydroSolutionType::Deflagration 
    ) {
        throw std::runtime_error(
            "The solution with v_+ = " + std::to_string(v_plus) + " and v_- = " 
            + std::to_string(v_minus) + " is not a deflagration but of type " 
            + to_string(classify_wave_wall_frame(v_plus, v_minus)) + "."
        );
    }
    // calculate fluid velocity in front / behind bubble in the fluid frame
    double xi_wall = v_minus;
    double v0 = mu(xi_wall, v_plus);
    if ( debug ) { 
        std::cout << "[DEBUG: compute_velocity_profile] Direct Deflagration:"
                  << " xi_w = " << xi_wall << ", v(xi_w) = " << v0 << "." 
                  << std::endl; 
    }
    
    // solve fluid equations
    int status; // integration return status
    xi.clear(); v.clear();
    // add first two elements
    xi.push_back(0.); xi.push_back(xi_wall);
    v.push_back(0.); v.push_back(0.);
    // set wall velocity index
    size_t idx_wall = xi.size(); // index of first element outside the bubble
    // add initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0);
    // integrate from wall to shock
    status = integrate_velocity_profile(
        xi, v, v_zero, step_size, true, solve_mode
    );
    if ( status != 1 ) {
        utils::check_failed(
            "[CHECKS: compute_velocity_profile_direct_deflagration] "
            "Expected integration return status 1 (encountered shock), but got "
            + std::to_string(status) + " (" 
            + (status==0 ? "reached v_stop" : "entered unphysical region") 
            + ") for v_+ = " + std::to_string(v_plus) + ", v_- = " 
            + std::to_string(v_minus) + "."
        );
        if ( xi.back() < cs ) {
            if ( debug ) {
                std::cout << "[DEBUG: "
                          << "compute_velocity_profile_direct_deflagration"
                          << "] encountered subsonic shock front at (" 
                          << xi.back() << ", " << v.back() 
                          << ") for direct deflagration with v_+ = " << v_plus
                          << " and v_- = " << v_minus << ". Appending (cs,0)." 
                          << std::endl;
            }
            xi.push_back(cs); v.push_back(0.);
        }
    }
    // add last two elements
    xi.push_back(xi.back()); xi.push_back(1.);
    v.push_back(0.); v.push_back(0.);
    
    // return index of wall velocity
    return idx_wall;
}

size_t compute_velocity_profile_direct_hybrid(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double xi_wall, 
    double step_size,
    ProfileSolveMode solve_mode
) {
    // check solution type
    double v_minus = cs;
    if ( 
        classify_wave_wall_frame(v_plus, v_minus) != HydroSolutionType::Hybrid 
    ) {
        throw std::runtime_error(
            "The solution with v_+ = " + std::to_string(v_plus) + " and v_- = " 
            + std::to_string(v_minus) + " is not a hybrid but of type " 
            + to_string(classify_wave_wall_frame(v_plus, v_minus)) + "."
        );
    }

    // calculate fluid velocity in front / behind bubble in the fluid frame
    double v0m = mu(xi_wall, v_minus);
    double v0p = mu(xi_wall, v_plus);
    if ( debug ) { 
        std::cout << "[DEBUG: compute_velocity_profile] Direct Hybrid: xi_w = " 
                  << xi_wall << ", v(xi_w^-) = " << v0m << ", v(xi_w^+) = " 
                  << v0p << "." << std::endl;
    }
    
    // solve fluid equations
    int status; // integration return status
    xi.clear(); v.clear();
    // initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0m);
    // integrate from wall to sound speed
    status = integrate_velocity_profile(
        xi, v, v_zero, step_size, false, solve_mode
    );
    if ( status != 0 ) {
        utils::check_failed(
            "[CHECKS: compute_velocity_profile_direct_hybrid] "
            "Expected integration return status 0 (reached v_stop), but got " 
            + std::to_string(status) + " (" + 
            (status==1 ? "encountered shock" : "entered unphysical region") 
            + ") for v_+ = " + std::to_string(v_plus) + ", xi_w = " 
            + std::to_string(xi_wall) + "."
        );
    }
    // add first two elements
    xi.push_back(cs); xi.push_back(0.);
    v.push_back(0.); v.push_back(0.);
    // revert order 
    std::reverse(xi.begin(), xi.end());
    std::reverse(v.begin(), v.end());
    // set wall velocity index
    size_t idx_wall = xi.size(); // index of first element outside the bubble
    // add initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0p);
    // integrate from wall to shock
    status = integrate_velocity_profile(
        xi, v, v_zero, step_size, true, solve_mode
    );
    if ( status != 1 ) {
        utils::check_failed(
            "[CHECKS: compute_velocity_profile_direct_hybrid] "
            "Expected integration return status 1 (encountered shock), but got "
            + std::to_string(status) + " (" 
            + (status==0 ? "reached v_stop" : "entered unphysical region") 
            + ") for v_+ = " + std::to_string(v_plus) + ", xi_w = " 
            + std::to_string(xi_wall) + "."
        );
    }
    // add last two elements
    xi.push_back(xi.back()); xi.push_back(1.);
    v.push_back(0.); v.push_back(0.);
    
    // return index of wall velocity
    return idx_wall;
}

size_t compute_velocity_profile_inverse_detonation(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double v_minus,
    double step_size,
    ProfileSolveMode solve_mode
) {
    // check solution type
    if ( 
        classify_wave_wall_frame(v_plus, v_minus) 
        != HydroSolutionType::InverseDetonation 
    ) {
        throw std::runtime_error(
            "The solution with v_+ = " + std::to_string(v_plus) + " and v_- = "
            + std::to_string(v_minus) + " is not an inverse detonation but of "
            "type " + to_string(classify_wave_wall_frame(v_plus, v_minus)) + "."
        );
    }
    // calculate fluid velocity in front / behind bubble in the fluid frame
    double xi_wall = v_minus;
    double v0 = mu(xi_wall, v_plus);
    if ( debug ) { 
        std::cout << "[DEBUG: compute_velocity_profile] Inverse Detonation:"
                  << " xi_w = " << xi_wall << " and v(xi_w) = " << v0 << "." 
                  << std::endl; 
    }
    
    // solve fluid equations
    int status; // integration return status
    xi.clear(); v.clear();
    // add first two elements
    xi.push_back(0.); xi.push_back(xi_wall);
    v.push_back(0.); v.push_back(0.);
    // set wall velocity index
    size_t idx_wall = xi.size(); // index of first element outside the bubble
    // add initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0);
    // integrate from wall to sound speed
    status = integrate_velocity_profile(
        xi, v, -v_zero, step_size, false, solve_mode
    );
    if ( status != 0 ) {
        utils::check_failed(
            "[CHECKS: compute_velocity_profile_inverse_detonation] "
            "Expected integration return status 0 (reached v_stop), but got " 
            + std::to_string(status) + " (" 
            + (status==1 ? "encountered shock" : "entered unphysical region") 
            + ") for v_+ = " + std::to_string(v_plus) + ", v_- = " 
            + std::to_string(v_minus) + "."
        );
    }
    // add last two elements
    xi.push_back(cs); xi.push_back(1.);
    v.push_back(0.); v.push_back(0.);
    
    // return index of wall velocity
    return idx_wall;
}

size_t compute_velocity_profile_inverse_deflagration(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double v_minus, 
    double step_size,
    ProfileSolveMode solve_mode
) {
    // check solution type
    if ( 
        classify_wave_wall_frame(v_plus, v_minus) != 
        HydroSolutionType::InverseDeflagration 
    ) {
        throw std::runtime_error(
            "The solution with v_+ = " + std::to_string(v_plus) + " and v_- = " 
            + std::to_string(v_minus) + " is not an inverse deflagration but of"
            " type " + to_string(classify_wave_wall_frame(v_plus, v_minus)) +"."
        );
    }

    // calculate fluid velocity in front / behind bubble in the fluid frame
    double xi_wall = v_plus;
    double v0 = mu(xi_wall, v_minus);
    if ( debug ) { 
        std::cout << "[DEBUG: compute_velocity_profile] Inverse Deflagration:"
                  << " xi_w = " << xi_wall << " and v(xi_w) = " << v0 << "." 
                  << std::endl; 
    }
    
    // solve fluid equations
    int status; // integration return status
    xi.clear(); v.clear();
    // add initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0);
    // integrate from wall to shock
    status = integrate_velocity_profile(
        xi, v, -v_zero, step_size, true, solve_mode
    );
    if ( status != 1 ) {
        utils::check_failed(
            "[CHECKS: compute_velocity_profile_inverse_deflagration] "
            "Expected integration return status 1 (encountered shock), but got "
            + std::to_string(status) + " (" 
            + (status==0 ? "reached v_stop" : "entered unphysical region") 
            + ") for v_+ = " + std::to_string(v_plus) + ", v_- = " 
            + std::to_string(v_minus) + "."
        );
        if ( xi.back() > cs ) {
            if ( debug ) {
                std::cout << "[DEBUG:"
                          << " compute_velocity_profile_inverse_deflagration]"
                          << " encountered supersonic shock front at (" 
                          << xi.back() << ", " << v.back() << ") for inverse"
                          << " deflagration with v_+ = " << v_plus 
                          << " and v_- = " << v_minus << ". Appending (cs,0)." 
                          << std::endl;
            }
            xi.push_back(cs); v.push_back(0.);
        }
    }
    // add first two elements
    xi.push_back(xi.back()); xi.push_back(0.);
    v.push_back(0.); v.push_back(0.);
    // revert order 
    std::reverse(xi.begin(), xi.end());
    std::reverse(v.begin(), v.end());
    // set wall velocity index
    size_t idx_wall = xi.size(); // index of first element outside the bubble
    // add last two elements
    xi.push_back(xi.back()); xi.push_back(1.);
    v.push_back(0.); v.push_back(0.);
    
    // return index of wall velocity
    return idx_wall;
}

size_t compute_velocity_profile_inverse_hybrid(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double xi_wall, 
    double v_minus, 
    double step_size,
    ProfileSolveMode solve_mode
) {
    // check solution type
    double v_plus = cs;
    if ( 
        classify_wave_wall_frame(v_plus, v_minus) 
        != HydroSolutionType::InverseHybrid 
    ) {
        throw std::runtime_error(
            "The solution with v_+ = " + std::to_string(v_plus) + " and v_- = "
            + std::to_string(v_minus) + " is not an inverse hybrid but of type "
            + to_string(classify_wave_wall_frame(v_plus, v_minus)) + "."
        );
    }
    // calculate fluid velocity in front / behind bubble in the fluid frame
    double v0m = mu(xi_wall, v_minus);
    double v0p = mu(xi_wall, v_plus);
    if ( debug ) { 
        std::cout << "[DEBUG: compute_velocity_profile] Inverse Hybrid: xi_w = "
                  << xi_wall << ", v(xi_w^-) = " << v0m << ", v(xi_w^+) = " 
                  << v0p << "." << std::endl;
    }

    // solve fluid equations
    int status; // integration return status
    xi.clear(); v.clear();
    // add initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0m);
    // integrate from wall to shock
    status = integrate_velocity_profile(
        xi, v, -v_zero, step_size, true, solve_mode
    );
    if ( status != 1 ) {
        utils::check_failed(
            "[CHECKS: compute_velocity_profile_inverse_hybrid] "
            "Expected integration return status 1 (encountered shock), but got "
            + std::to_string(status) + " (" 
            + (status==0 ? "reached v_stop" : "entered unphysical region") 
            + ") for xi_w = " + std::to_string(xi_wall) + ", v_- = " 
            + std::to_string(v_minus) + "."
        );
    }
    // add first two elements
    xi.push_back(xi.back()); xi.push_back(0.);
    v.push_back(0.); v.push_back(0.);
    // revert order 
    std::reverse(xi.begin(), xi.end());
    std::reverse(v.begin(), v.end());
    // set wall velocity index
    size_t idx_wall = xi.size(); // index of first element outside the bubble
    // add initial conditions
    xi.push_back(xi_wall);
    v.push_back(v0p);
    // integrate from wall to sound speed
    status = integrate_velocity_profile(
        xi, v, -v_zero, step_size, false, solve_mode
    );
    if ( status != 0 ) {
        utils::check_failed(
            "[CHECKS: compute_velocity_profile_inverse_hybrid] "
            "Expected integration return status 0 (reached v_stop), but got " 
            + std::to_string(status) + " (" 
            + (status==1 ? "encountered shock" : "entered unphysical region") 
            + ") for xi_w = " + std::to_string(xi_wall) + ", v_- = " 
            + std::to_string(v_minus) + "."
        );
    }
    // add last two elements
    xi.push_back(cs); xi.push_back(1.);
    v.push_back(0.); v.push_back(0.);
    
    // return index of wall velocity
    return idx_wall;
}

size_t get_shock_index(HydroSolutionType sol_type, size_t size) {
    switch (sol_type) {
        case HydroSolutionType::Deflagration:
        case HydroSolutionType::Hybrid:
            return size - 2;
        case HydroSolutionType::InverseDeflagration:
        case HydroSolutionType::InverseHybrid:
            return 2;
        default:
            break;
    }
    throw std::runtime_error(
        "Solution type " + to_string(sol_type) + " does not have a shock front."
    );
}

std::pair<double, double> get_shock_front(
    HydroSolutionType sol_type, 
    const std::vector<double>& xi, 
    const std::vector<double>& v
) {
    check_profiles(xi, v);
    size_t idx = get_shock_index(sol_type, xi.size());
    if ( 
        sol_type == HydroSolutionType::Deflagration 
        || sol_type == HydroSolutionType::Hybrid 
    ) { --idx; } // for direct transition, return v_-, else v_+
    return {xi[idx], v[idx]};
}

size_t compute_velocity_profile(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double xi_wall, 
    double alpha_N,
    double step_size,
    ProfileSolveMode solve_mode
) {
    HydroSolutionType sol_type = classify_wave_fluid_frame(
        xi_wall, alpha_N, step_size, solve_mode
    );
    double alpha_plus = find_alpha_plus(
        sol_type, alpha_N, xi_wall, step_size, solve_mode
    );
    return compute_velocity_profile(
        sol_type, xi, v, xi_wall, alpha_plus, step_size, solve_mode
    );
}

size_t compute_velocity_profile(
    HydroSolutionType sol_type, 
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double xi_wall, 
    double alpha_p, 
    double step_size,
    ProfileSolveMode solve_mode
) {
    
    // solve matching conditions at bubble wall
    auto [v_plus, v_minus] = compute_fluid_velocities(
        xi_wall, alpha_p, sol_type
    );

    size_t idx_wall; // index of first element outside the bubble
    switch (sol_type) {
        case HydroSolutionType::Deflagration:
            idx_wall = compute_velocity_profile_direct_deflagration(
                xi, v, v_plus, v_minus, step_size, solve_mode
            );
            break;
        case HydroSolutionType::Detonation:
            idx_wall = compute_velocity_profile_direct_detonation(
                xi, v, v_plus, v_minus, step_size, solve_mode
            );
            break;
        case HydroSolutionType::Hybrid:
            idx_wall = compute_velocity_profile_direct_hybrid(
                xi, v, v_plus, xi_wall, step_size, solve_mode
            );
            break;
        case HydroSolutionType::InverseDeflagration:
            idx_wall = compute_velocity_profile_inverse_deflagration(
                xi, v, v_plus, v_minus, step_size, solve_mode
            );
            break;
        case HydroSolutionType::InverseDetonation:
            idx_wall = compute_velocity_profile_inverse_detonation(
                xi, v, v_plus, v_minus, step_size, solve_mode
            );
            break;
        case HydroSolutionType::InverseHybrid:
            idx_wall = compute_velocity_profile_inverse_hybrid(
                xi, v, xi_wall, v_minus, step_size, solve_mode
            );
            break;
        default:
            throw std::runtime_error(
                "Unknown solution type: " + to_string(sol_type) + "."
            );
    }

    if ( !are_equal(xi_wall, xi[idx_wall]) ) {
        utils::check_failed(
            "[CHECKS: compute_velocity_profile] "
            "The wall velocity of the calculated velocity spectrum does not"
            " match the input value: " + std::to_string(xi[idx_wall]) + " vs. " 
            + std::to_string(xi_wall) + "."
        );
    }

    // check shock front
    double v_shock = v[2];      // for inverse deflagration/hybrid
    double dv = v_shock - v[3]; // for inverse deflagration/hybrid
    switch (sol_type) {
        case HydroSolutionType::Deflagration:
        case HydroSolutionType::Hybrid:
            v_shock = v[v.size()-3];
            dv = v[v.size()-4] - v_shock;
        case HydroSolutionType::InverseDeflagration:
        case HydroSolutionType::InverseHybrid:
            if ( std::abs(v_shock) < std::abs(dv) ) {
                utils::check_failed(
                    "[CHECKS: compute_velocity_profile] v_shock " 
                    + std::to_string(v_shock) 
                    + " is too small or zero. Increase precision!"
                );
            }
            break;
        case HydroSolutionType::Detonation:
        case HydroSolutionType::InverseDetonation:
        default:
            break;
    }

    return idx_wall;
}

// TODO: can be split into constant and linear part. Constant part needs to be 
//       evaluated at start/end and jumps only.
double integrate_enthalpy_segment(double xi1, double xi2, double v1, double v2)
{
    if ( xi1 == xi2 ) { return 0.0; } // jump position
    if ( v1 == v2 ) { return 0.0; } // integral vanishes if v is constant
    // approximate xi(v) as linear function: xi(v) = xi0 + dxi * v
    double xi0 = (xi1 * v2 - xi2 * v1) / (v2 - v1);
    double dxi = (xi2 - xi1) / (v2 - v1);
    // partial fraction decomposition:
    // gamma^2 * mu = (xi0-(1-dxi)*v) / (1-v^2) / (1-xi0*v-dxi*v^2) 
    //              = - v/(1-v^2) - (dxi*v + xi0)/(dxi*v^2+xi0*v-1))
    // first integral:
    double result = 0.5 * (std::log(1.0 - v2*v2) - std::log(1.0 - v1*v1));
    // second integral:
    result += 0.5 * (std::log(1.0 - xi1*v1) - std::log(1.0 - xi2*v2));
    double delta = 4.0 * dxi + xi0*xi0;
    if ( delta > 0 ) { 
        result += 0.5 * xi0 / std::sqrt(delta) * std::log(
            (2.0 * dxi * v1 + xi0 - std::sqrt(delta)) 
            / (2.0 * dxi * v1 + xi0 + std::sqrt(delta)) 
            * (2.0 * dxi * v2 + xi0 + std::sqrt(delta))
            / (2.0 * dxi * v2 + xi0 - std::sqrt(delta))
        );
    } else {
        result += xi0 / std::sqrt(-delta) * (
            std::atan((2.0 * dxi * v1 + xi0) / std::sqrt(-delta))
            - std::atan((2.0 * dxi * v2 + xi0) / std::sqrt(-delta))
        );
    }
    return (1.0 + 1.0/cs2) * result;
}

std::vector<double> compute_enthalpy_density(
    const std::vector<double>& xi, 
    const std::vector<double>& v, 
    size_t idx_wall, 
    size_t idx_shock
) {
    check_profiles(xi, v);
    const size_t N = xi.size();
    std::vector<double> w(N);
    w[N-1] = 1.0;
    for ( size_t i=xi.size()-1; i>0; --i ) {
        // note: if idx_shock=-1, its actual value is PTRDIFF_MAX-1
        //       (since size_t is unsigned) but this is fine here
        if ( i == idx_wall || i == idx_shock ) { 
            if ( xi[i] != xi[i-1] ) { 
                throw std::runtime_error(
                    std::to_string(i) + " is not a jump position: No jump"
                    " between (xi,v) = (" + std::to_string(xi[i-1]) + ", " 
                    + std::to_string(v[i-1]) + ") and (xi, v) = ()"
                    + std::to_string(xi[i]) + ", " + std::to_string(v[i]) + "."
                );
            }
            // jump in enthalpy density: go to rest frame of the front and match
            double vp = mu(xi[i], v[i]);
            double vm = mu(xi[i], v[i-1]);
            if ( !check_velocity(vp) || !check_velocity(vm) ) {
                throw std::runtime_error(
                    "Encountered invalid velocities after boosting back to wall"
                    " frame: v_+ = " + std::to_string(vp) + " and v_- = " 
                    + std::to_string(vm) + "."
                );
            }
            if ( debug ) { 
                std::cout << "[DEBUG: compute_enthalpy_density] reconstructed"
                          << " (v_+, v-)  = (" << vp << ", " << vm << ") at "
                          << (i==idx_wall ? "bubble" : "shock") << " front." 
                          << std::endl;
            }
            w[i-1] = w[i] * vp / vm * (1.0 - vm * vm) / (1.0 - vp * vp);
        } else {
            w[i-1] = w[i] * std::exp(integrate_enthalpy_segment(
                xi[i], xi[i-1], v[i], v[i-1]
            ));
        }
    }
    return w;
}

std::vector<double> compute_energy_density_profile(
    const std::vector<double>& w, size_t idx_wall, double epsilon
) {
    std::vector<double> e(w.size());
    for ( size_t i=0; i<w.size(); ++i )  { 
        e[i] = 0.75 * w[i] + (i < idx_wall ? 0.0 : epsilon); 
    }
    return e;
}

std::vector<double> compute_pressure_density_profile(
    const std::vector<double>& w, size_t idx_wall, double epsilon
) {
    std::vector<double> p(w.size());
    for ( size_t i=0; i<w.size(); ++i )  { 
        p[i] = 0.25 * w[i] - (i < idx_wall ? 0.0 : epsilon); 
    }
    return p;
}

double compute_epsilon(double w, double alpha) { return 0.75 * alpha * w; }

double compute_alpha_N(
    HydroSolutionType sol_type, 
    double alpha_plus, 
    double xi_wall, 
    double step_size,
    ProfileSolveMode solve_mode
) {
    auto [v_plus, v_minus] = compute_fluid_velocities(
        xi_wall, alpha_plus, sol_type
    ); 
    std::vector<double> xi{xi_wall}, v{mu(xi_wall, v_plus)};
    double xish, mush, alpha_N;
    double matching_factor = 1.0;
    switch (sol_type) {
        // for detonation and inverse deflagration:
        //     alpha_N and alpha_+ are identical
        case HydroSolutionType::Detonation:
        case HydroSolutionType::InverseDeflagration:
            return alpha_plus;
        // for deflagration and hybrids: 
        //     evolve w from bubble wall to shock front
        case HydroSolutionType::Deflagration:
        case HydroSolutionType::Hybrid:
            integrate_velocity_profile(
                xi, v, v_zero, step_size, true, solve_mode
            );
            xish = xi.back(); // vsh+
            mush = mu(xish, v.back()); // vsh-
            if ( !check_velocity(mush) ) {
                throw std::runtime_error(
                    "Encountered invalid velocity after boosting back to shock"
                    " frame: mu(xi_sh, v_sh) = " + std::to_string(mush) + "."
                );
            }
            matching_factor =  (
                xish / mush * (1.0 - mush * mush) / (1.0 - xish*xish)
            ); // matching at shock front
            break;
        // for inverse detonation and inverse hybrids: evolve w from bubble wall
        // to speed-of-sound
        case HydroSolutionType::InverseHybrid:
            v[0] = mu(xi_wall, cs);
            // no break here!
        case HydroSolutionType::InverseDetonation:
            integrate_velocity_profile(
                xi, v, -v_zero, step_size, false, solve_mode
            );
            break;
    }
    double log_w_N_over_w_plus = 0.0;
    for ( size_t i=0; i<xi.size()-1; ++i ) { 
        log_w_N_over_w_plus += integrate_enthalpy_segment(
            xi[i+1], xi[i], v[i+1], v[i]
        ); 
    }
    alpha_N = alpha_plus * matching_factor * std::exp(log_w_N_over_w_plus);
    return alpha_N;
}

double find_alpha_plus(
    double alpha_N, 
    double xi_wall, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    HydroSolutionType sol_type = classify_wave_fluid_frame(
        xi_wall, alpha_N, step_size, solve_mode
    ); 
    return find_alpha_plus(sol_type, alpha_N, xi_wall, step_size, solve_mode);
}

double find_alpha_plus(
    HydroSolutionType sol_type,
    double alpha_N, 
    double xi_wall, 
    double step_size, 
    ProfileSolveMode solve_mode
) {    
    // bracket alpha_+ based on solution type
    // (initialize to alpha_N to avoid compiler warnings)
    double alpha_min = alpha_N, alpha_max = alpha_N;
    switch ( sol_type ) {
        // note: for direct detonation and inverse deflagration, alpha_N and
        //       alpha_+ are the same
        case HydroSolutionType::Detonation:
        case HydroSolutionType::InverseDeflagration:
            if ( 
                alpha_N < alpha_min_inverse_deflagration(xi_wall, false) || 
                alpha_N > alpha_max_detonation(xi_wall, false) ||
                xi_wall < cs
            ) { throw std::runtime_error(
                "Invalid combination of wall velocity and transition strength:"
                " xi_w = " + std::to_string(xi_wall) + ", alpha_N = " 
                + std::to_string(alpha_N) + " for " + to_string(sol_type) 
                + " solutions."
            ); }
            return alpha_N; 
        case HydroSolutionType::Deflagration:
            alpha_min = alpha_abs_tol;
            // 1.0 / 3.0
            alpha_max = alpha_max_deflagration(xi_wall, false) - alpha_abs_tol; 
            break;
        case HydroSolutionType::Hybrid:
            // hybrids cannot change to detonations!
            alpha_min = alpha_min_hybrid(xi_wall, false) + alpha_abs_tol; 
            // 1.0 / 3.0
            alpha_max = alpha_max_hybrid(xi_wall, false) - alpha_abs_tol;
            break;
        case HydroSolutionType::InverseDetonation:
            alpha_min = (
                alpha_min_inverse_detonation(xi_wall, false) + alpha_abs_tol
            );
            alpha_max = -alpha_abs_tol;
            break;
        case HydroSolutionType::InverseHybrid:
            alpha_min = (
                alpha_min_inverse_hybrid(xi_wall, false) + alpha_abs_tol
            );
            alpha_max = (
                alpha_min_inverse_detonation(xi_wall, false) - alpha_abs_tol
            );
            break;
    }
    alpha_max = std::min(alpha_max, alpha_N);
    
    if ( debug ) {
        std::cout << "[DEBUG: find_alpha_plus] " << sol_type << ": bracketed"
                  << " alpha in [" << alpha_min << ", " << alpha_max << "] for"
                  << " (xi_w, alpha_N) = " << xi_wall << ", " << alpha_N << ")."
                  << std::endl;
    }

    // find root of the relative difference
    auto alpha_rel_difference = [
        sol_type, alpha_N, xi_wall, step_size, solve_mode
    ](double alpha_plus) { 
        return compute_alpha_N(
            sol_type, alpha_plus, xi_wall, step_size, solve_mode
        ) / alpha_N - 1.0; 
    };
    // use combination of absolute and relative tolerance
    auto tolerance = [
        rel_tol = alpha_rel_tol, abs_tol = alpha_abs_tol
    ](double a, double b) {
        return std::abs(b - a) <= std::max(
            abs_tol, 
            rel_tol*std::max(std::abs(a), std::abs(b))
        );
    };

    // check that alpha_N really changes sign in the bracketing interval
    if ( 
        alpha_rel_difference(alpha_max) * alpha_rel_difference(alpha_min) > 0. 
    ) {
        throw std::runtime_error(
            "Cannot find alpha_+: No sign change in initial bracketing range ("
            + std::to_string(alpha_min) + ", " + std::to_string(alpha_max) 
            + "): alpha_N = (" + std::to_string(compute_alpha_N(
                sol_type, alpha_min, xi_wall, step_size, solve_mode
            )) + ", " + std::to_string(compute_alpha_N(
                sol_type, alpha_max, xi_wall, step_size, solve_mode
            )) + ")."
        );
    }

    // find root
    std::uintmax_t N_iter = 1000;
    auto result = boost::math::tools::toms748_solve(
        alpha_rel_difference, alpha_min, alpha_max, tolerance, N_iter
    );
    double alpha_plus = .5*(result.first+result.second);

    if ( debug ) {
        double delta_alpha_plus = alpha_rel_difference(alpha_plus);
        std::cout << "[DEBUG: find_alpha_plus] Found alpha_+ = " << alpha_plus 
                  << ". " << "Relative difference in alpha_N: " 
                  << delta_alpha_plus << "." << std::endl;
    }

    return alpha_plus;
}

double compute_efficiency_factor(
    const std::vector<double>& xi,  
    const std::vector<double>& v, 
    const std::vector<double>& w,
    size_t idx_wall, 
    double epsilon
) {
    size_t size = xi.size();
    check_profiles(xi, v); check_profiles(xi, w);
    std::vector<double> integrand(size);
    double xi_wall = xi[idx_wall];
    double prefactor = (
        epsilon > 0.
        ?
        3./epsilon / xi_wall/xi_wall/xi_wall
        :
        4./w.back() / std::max(xi_wall*xi_wall*xi_wall, cs*cs2)
    );
    #pragma omp simd
    for ( size_t i=0; i<size; ++i ) {
        double temp = v[i] * gamma(v[i]) * xi[i];
        integrand[i] = w[i] * temp*temp;
    }
    return prefactor * utils::integrate_trapez(xi, integrand);
}

double compute_dofs_ratio_LTE(
    double xi_wall, 
    double alpha_N, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    HydroSolutionType sol_type = classify_wave_fluid_frame(
        xi_wall, alpha_N, step_size, solve_mode
    );
    double alpha_plus = find_alpha_plus(
        sol_type, alpha_N, xi_wall, step_size, solve_mode
    );
    auto [v_plus, v_minus] = compute_fluid_velocities(
        xi_wall, alpha_plus, sol_type
    );
    if ( sol_type == HydroSolutionType::Hybrid ) { v_minus = cs; }
    else if ( sol_type == HydroSolutionType::InverseHybrid ) { v_plus = cs; }
    return v_plus*(1.0 - v_plus*v_plus) / v_minus/(1.0 - v_minus*v_minus);
}

double find_LTE_wall_velocity(
    double alpha_N, 
    double Psi, 
    double xi_min, 
    double xi_max, 
    double step_size, 
    ProfileSolveMode solve_mode
) {
    if ( xi_min >= xi_max ) { throw std::runtime_error(
        "Cannot find LTE wall velocity: xi_min and xi_max are not properly"
        " ordered: xi_min = " + std::to_string(xi_min) + " > xi_max = " 
        + std::to_string(xi_max) + "."
    ); }
    if ( 
        alpha_N > 0 
        && xi_min < direct_jouguet_velocity(alpha_N) 
        && xi_max > direct_jouguet_velocity(alpha_N) 
    ) {
        throw std::runtime_error(
            "Cannot find LTE wall velocity: The bracketing inverval (" 
            + std::to_string(xi_min) + ", " + std::to_string(xi_max) 
            + ") crosses the Jouguet velocity " 
            + std::to_string(direct_jouguet_velocity(alpha_N)) + "."
        );
    } else if ( alpha_N < 0 && xi_min < cs && xi_max > cs ) {
        throw std::runtime_error(
            "Cannot find LTE wall velocity: The bracketing inverval (" 
            + std::to_string(xi_min) + ", " + std::to_string(xi_max) 
            + ") crosses the speed of sound."
        );
    }
    if ( debug ) {
        std::cout << "[DEBUG: find_LTE_wall_velocity] Bracketing LTE wall"
                  << " velocity in [" << xi_min << ", " << xi_max << "] for"
                  << " (alpha_N, Psi) = (" << alpha_N << ", " << Psi << ")."
                  << std::endl;
    }
    // find root of the relative difference
    auto psi_rel_difference = [alpha_N, Psi, step_size, solve_mode](
        double xi_wall
    ) { 
        return compute_dofs_ratio_LTE(
            xi_wall, alpha_N, step_size, solve_mode
        ) / Psi - 1.0; 
    };
    // use combination of absolute and relative tolerance
    auto tolerance = [rel_tol = xi_wall_rel_tol, abs_tol = xi_wall_abs_tol](
        double a, double b
    ) {
        return std::abs(b - a) <= std::max(
            abs_tol, 
            rel_tol*std::max(std::abs(a), std::abs(b))
        );
    };
    // check that Psi really changes sign in the bracketing interval
    if ( psi_rel_difference(xi_min) * psi_rel_difference(xi_max) > 0. ) {
        throw std::runtime_error(
            "Cannot find LTE wall velocity: No sign change in initial"
            " bracketing range (" + std::to_string(xi_min) + ", " 
            + std::to_string(xi_max) + "): Delta_Psi/Psi = ("
            + std::to_string(psi_rel_difference(xi_min)) + ", " 
            + std::to_string(psi_rel_difference(xi_max)) + ")."
        );
    }

    // find root
    std::uintmax_t N_iter = 1000;
    auto result = boost::math::tools::toms748_solve(
        psi_rel_difference, xi_min, xi_max, tolerance, N_iter
    );
    double xi_wall = .5*(result.first+result.second);

    if ( debug ) {
        double delta_Psi = psi_rel_difference(xi_wall);
        std::cout << "[DEBUG: find_LTE_wall_velocity] Found xi_wall = " 
                  << xi_wall << ". " << "Relative difference in Psi: " 
                  << delta_Psi << "." << std::endl;
    }

    return xi_wall;
}

bool check_profiles(
    const std::vector<double>& xi, 
    const std::vector<double>& profile, 
    bool throw_error
) {
    size_t size = xi.size();
    if ( size != profile.size() || size < 6) {
        std::string message = 
            "The self-similar coordinate and the profile vectors must have the"
            " same size and at least size 6. Got sizes " + std::to_string(size) 
            + " and " + std::to_string(profile.size()) + ".";

        if ( throw_error ) { throw std::runtime_error(message); }
        if ( debug ) { 
            std::cerr << "[DEBUG: check_profiles] " << message << std::endl; 
        }
        return false;
    }
    return true;
}

} // namespace profile_calculator
