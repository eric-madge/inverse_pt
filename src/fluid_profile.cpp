/*
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

#include "inverse_pt/fluid_profile.hpp"

#include "inverse_pt/constants.hpp"
#include "inverse_pt/profile_calculator.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace inverse_pt {

FluidProfile::FluidProfile(
    const std::vector<double>& xi_, 
    const std::vector<double>& v_, 
    const std::vector<double>& w_, 
    const std::vector<double>& e_, 
    const std::vector<double>& p_,
    size_t idx_w, 
    size_t idx_sh
) : 
    xi(xi_), 
    v(v_), 
    w(w_), 
    e(e_), 
    p(p_), 
    idx_wall(idx_w), 
    idx_shock(idx_sh), 
    xi_wall(xi_[idx_w]) 
{ 
    if ( verbose ) { 
        std::cout << "Initializing fluid profile from data. " << std::endl; 
    }
    // check profiles
    check_profiles(); 
    if ( verbose ) { std::cout << "Profiles are okay." << std::endl; }
    // classify solution
    double v_plus = profile_calculator::mu(xi_wall, v[idx_wall]);
    double v_minus = profile_calculator::mu(xi_wall, v[idx_wall-1]);
    sol_type = profile_calculator::classify_wave_wall_frame(v_plus, v_minus);
    double alpha_plus = profile_calculator::compute_alpha_plus(v_plus, v_minus);
    double alpha_N = alpha_plus * w[idx_wall] / w.back();
    alpha = alpha_N;
    if ( verbose ) {
        std::cout << "Classified profile as " << sol_type << " solution with"
                  << " v_+ = " << v_plus << " and v_- = " << v_minus << ","
                  << " alpha_+ = " << alpha_plus << ".\nThis corresponds to"
                  << " v_w = " << xi_wall << " and alpha_N = " << alpha_N << "."
                  << std::endl; 
    }
    if ( ( idx_shock == constants::NO_INDEX ) != ( 
        sol_type == profile_calculator::HydroSolutionType::Detonation || 
        sol_type == profile_calculator::HydroSolutionType::InverseDetonation 
    ) ) { throw std::runtime_error(
        "The solution is of type " + to_string(sol_type) + " but " + ( 
            sol_type == profile_calculator::HydroSolutionType::Detonation || 
            sol_type == profile_calculator::HydroSolutionType::InverseDetonation 
            ? " shock index " + std::to_string(idx_shock) : " no shock index "
        ) + " was provided."
    ); } 
}

FluidProfile::FluidProfile(
    double v_w, 
    double alpha_N, 
    double step_size, 
    profile_calculator::ProfileSolveMode solve_mode
) : xi_wall(v_w), alpha(alpha_N) {
    if ( verbose ) {
        std::cout << "Calculating fluid profile for xi_w = " << v_w 
                  << " and alpha_N = " << alpha_N << "." << std::endl;
    }

    sol_type = profile_calculator::classify_wave_fluid_frame(
        xi_wall, alpha_N, step_size, solve_mode
    );
    double alpha_plus = profile_calculator::find_alpha_plus(
        sol_type, alpha_N, xi_wall, step_size, solve_mode
    );
    auto [vplus, vminus] = profile_calculator::compute_fluid_velocities(
        xi_wall, alpha_plus, sol_type
    );

    if ( verbose ) {
        std::cout << "Found " << sol_type << " solution with alpha_+ = " 
                  << alpha_plus << ".\n" << "Fluid velocities: v+ = " << vplus 
                  << ", v_- = " << vminus << ".\nCalculating velocity profile."
                  << std::endl;
    }

    idx_wall = profile_calculator::compute_velocity_profile(
        sol_type, xi,v, xi_wall, alpha_plus, step_size, solve_mode
    );
    idx_shock = constants::NO_INDEX;
    if ( 
        sol_type != profile_calculator::HydroSolutionType::Detonation 
        && sol_type != profile_calculator::HydroSolutionType::InverseDetonation
    ) {
        idx_shock = profile_calculator::get_shock_index(sol_type, xi.size());
    }
    if ( verbose ) {
        std::cout << "Velocity profile calculated with " << xi.size() 
                  << " points.\n  Velocity jumps from " << v[idx_wall-1] 
                  << " to " << v[idx_wall] << " at the bubble wall, xi = " 
                  << xi[idx_wall] << " (index: " << idx_wall << ")." ;
        if ( idx_shock != constants::NO_INDEX ) { 
            std::cout << "\n  Velocity jumps from " << v[idx_shock-1] << " to " 
                      << v[idx_shock] << " at the shock position, xi = " 
                      << xi[idx_shock] << " (index: " << idx_shock << ").";
        }
        std::cout << std::endl;
    }
    w = profile_calculator::compute_enthalpy_density(
        xi, v, idx_wall, idx_shock
    );
    double wN = w.back();
    if (verbose) { 
        std::cout << "Enthalpy density calculated. wN: " << wN << std::endl; 
    }

    double epsilon = profile_calculator::compute_epsilon(wN, alpha_N);
    if ( verbose ) { 
        std::cout << "Vacuum energy: epsilon = " << epsilon << std::endl; 
    }
    e = profile_calculator::compute_energy_density_profile(
        w, idx_wall, epsilon
    );
    p = profile_calculator::compute_pressure_density_profile(
        w, idx_wall, epsilon
    );

    if ( verbose ) { std::cout << "All profiles calculated." << std::endl; }
}

const std::vector<double>& FluidProfile::get_xi_coordinate() const 
    { return xi; }
const std::vector<double>& FluidProfile::get_fluid_velocity() const 
    { return v; }
const std::vector<double>& FluidProfile::get_enthalpy_density() const 
    { return w; }
const std::vector<double>& FluidProfile::get_energy_density() const 
    { return e; }
const std::vector<double>& FluidProfile::get_pressure_density() const 
    { return p; }
size_t FluidProfile::get_wall_index() const { return idx_wall; }
size_t FluidProfile::get_shock_index() const { return idx_shock; }
double FluidProfile::get_wall_velocity() const { return xi_wall; }
double FluidProfile::get_transition_strength() const { return alpha; }
profile_calculator::HydroSolutionType FluidProfile::get_solution_type() const
    { return sol_type; }
bool FluidProfile::has_shock() const 
    { return idx_shock != constants::NO_INDEX; }
double FluidProfile::get_mean_enthalpy_density() const { return w.back(); }
double FluidProfile::get_adiabatic_index() const { return 4.0/3.0; }

double FluidProfile::get_efficiency_factor() const {
    double eps = profile_calculator::compute_epsilon(w.back(), alpha);
    return profile_calculator::compute_efficiency_factor(
        xi, v, w, idx_wall, eps
    );
}

void FluidProfile::set_verbose(bool b) { verbose = b; };
bool FluidProfile::is_verbose() { return verbose; };

void FluidProfile::check_profiles() const {
    profile_calculator::check_profiles(xi, v);
    profile_calculator::check_profiles(xi, w);
    profile_calculator::check_profiles(xi, e);
    profile_calculator::check_profiles(xi, p);
    if ( idx_wall < 2 || idx_wall > xi.size() - 2 ) { throw std::runtime_error(
        "Invalid value for the wall velocity index: " + std::to_string(idx_wall)
        + ". The index must be between 2 and " + std::to_string(xi.size()-2)
        + " (for size " + std::to_string(xi.size()) + ")." 
    );}
    if ( ! profile_calculator::are_equal(xi[idx_wall-1], xi[idx_wall]) ) { 
        throw std::runtime_error(
            "Wall index " + std::to_string(idx_wall) + " does not correspond to"
            " the bubble wall. " + "The self-similar coordinate changes from " 
            + std::to_string(xi[idx_wall-1]) + " to " 
            + std::to_string(xi[idx_wall]) + "."
        );
    }
    if ( idx_shock != constants::NO_INDEX ) {
        if ( idx_shock < 2 || idx_shock > xi.size() - 2 ) { 
            throw std::runtime_error(
                "Invalid value for the shock front index: " 
                + std::to_string(idx_shock) + ". The index must be between 2"
                " and " + std::to_string(xi.size()-2) + " (for size " 
                + std::to_string(xi.size()) + ")." 
            );
        }
        if ( ! profile_calculator::are_equal(xi[idx_shock-1], xi[idx_shock]) ) {
            throw std::runtime_error(
                "Shock index " + std::to_string(idx_shock) + " does not"
                " correspond to the shock front. The self-similar coordinate"
                " changes from " + std::to_string(xi[idx_shock-1]) + " to " + 
                std::to_string(xi[idx_shock]) + "."
            );
        }  
    }
}

} // namespace inverse_pt