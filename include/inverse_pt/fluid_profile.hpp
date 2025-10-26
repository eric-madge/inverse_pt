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

#include "inverse_pt/constants.hpp"
#include "inverse_pt/profile_calculator.hpp"

#include <vector>

namespace inverse_pt {

/**
 * @brief Class representing the fluid profile of a bubble in a cosmological
 *     phase transition.
 */
class FluidProfile {
public:

    /// @name Constructors
    ///@{

    /**
     * @brief Constructor for the \c FluidProfile from precomputed profiles.
     * @param xi_ Self-similar coordinate \f$\xi = r/t\f$
     * @param v_ Fluid velocity profile \f$v(\xi)\f$
     * @param w_ Enthalpy density profile \f$w(\xi)\f$
     * @param e_ Energy density profile \f$e(\xi)\f$
     * @param p_ Pressure density profile \f$p(\xi)\f$
     * @param idx_w Index of the first element outside the bubble
     * @param idx_sh Index of the first element outside the shock front
     *     (use \c constants::NO_INDEX if there is no shock front)
     */
    FluidProfile(
        const std::vector<double>& xi_, 
        const std::vector<double>& v_, 
        const std::vector<double>& w_, 
        const std::vector<double>& e_, 
        const std::vector<double>& p_, 
        size_t idx_w, 
        size_t idx_sh = constants::NO_INDEX
    );

    /** 
     * @brief Constructor for the \c FluidProfile from the wall velocity
     *     \f$\xi_w\f$ and strength parameter \f$\alpha_N\f$.
     * @param v_w wall veclocity \f$\xi_w\f$
     * @param alpha_N strength parameter \f$\alpha_N\f$
     *     (far in front of the wall)
     * @param step_size Step size for the fluid velocity (or coordinate) when
     *     solving the differential equation (default: 1e-3)
     * @param solve_mode \c ProfileSolveMode switch determining which
     *     differential equation for the velocity profile is solved (default:
     *     \c profile_calculator::ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
     */
    FluidProfile(
        double v_w, 
        double alpha_N, 
        double step_size=1e-3,
        profile_calculator::ProfileSolveMode solve_mode 
            = profile_calculator::ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
    );

    ///@}

    /** 
     * @brief Get the self-similar coordinate \f$\xi=r/t\f$ of the profile.
     * @return Reference to the vector containing the coordinate
     */
    const std::vector<double>& get_xi_coordinate() const;

    /** 
     * @brief Get the fluid velocity profile \f$v(\xi)\f$.
     * @return Reference to the vector containing the velocity
     */
    const std::vector<double>& get_fluid_velocity() const;

    /** 
     * @brief Get the enthalpy density profile \f$w(\xi)\f$.
     * @return Reference to the vector containing the enthalpy
     */
    const std::vector<double>& get_enthalpy_density() const;

    /** 
     * @brief Get the energy density profile \f$e(\xi)\f$.
     * @return Reference to the vector containing the energy
     */
    const std::vector<double>& get_energy_density() const;

    /** 
     * @brief Get the pressure density profile \f$p(\xi)\f$.
     * @return Reference to the vector containing the pressure
     */
    const std::vector<double>& get_pressure_density() const;

    /**
     * @brief Get the index of the element in \f$\xi\f$ corresponding to the
     *     bubble wall position.
     * @return Index of the bubble wall position (first element outside the
     *     wall)
     */
    size_t get_wall_index() const;

    /**
     * @brief Get the index of the element in \f$\xi\f$ corresponding to the
     *     shock front position.
     * @return Index of the shock front position (first element outside the
     *     shock)
     */
    size_t get_shock_index() const;

    /**
     * @brief Get the wall velocity \f$\xi_w\f$.
     * @return Value of the wall velocity \f$\xi_w\f$
     */
    double get_wall_velocity() const;

    /**
     * @brief Get the phase transition strength \f$\alpha_N\f$.
     * @return Value of the transition strength \f$\alpha_N\f$
     */
    double get_transition_strength() const;

    /**
     * @brief Get the type of the hydrodynamic solution.
     * @return Hydrodynamic solution type
     */
    profile_calculator::HydroSolutionType get_solution_type() const;

    /**
     * @brief Check whether the profile has a shock front. 
     * @return \c true if the profile has a shock front, \c false otherwise
     */
    bool has_shock() const;

    /**
     * @brief Get the mean enthalpy density \f$\bar{\omega}\f$ far in front of
     *     the wall.
     * @return Mean enthalpy density \f$\bar{w}\f$
     */
    double get_mean_enthalpy_density() const;

    /**
     * @brief Get the mean adiabatic index \f$\Gamma=\frac{\bar{w}}{\bar{e}}\f$
     *     in the stable phase.
     * @return Mean adiabatic index \f$\Gamma\f$
     */
    double get_adiabatic_index() const;

    /** 
     * @brief Get the efficiency factor \f$\kappa_\mathrm{sw}\f$.
     * @return Efficiency factor \f$\kappa_\mathrm{sw}\f$
     */
    double get_efficiency_factor() const;

    /**
     * @brief Set verbosity.
     * @param b Boolean value to which the verbosity flag is set
     */
    static void set_verbose(bool b);

    /**
     * @brief Get verbosity.
     * @return Boolean value to which the verbosity flag is set
     */
    static bool is_verbose();

private:
    inline static bool verbose = false; ///< verbosity flag
    std::vector<double> xi; ///< self-similar coordinate \f$\xi=r/t\f$
    std::vector<double> v;  ///< fluid velocity profile \f$v(\xi)\f$
    std::vector<double> w;  ///< enthalpy density profile \f$w(\xi)\f$
    std::vector<double> e;  ///< energy density profile \f$e(\xi)\f$
    std::vector<double> p;  ///< pressure density profile \f$p(\xi)\f$
    size_t idx_wall;        ///< index of the first element outside the wall
    size_t idx_shock;       ///< index of the first element outside the shock
    double xi_wall;         ///< wall velocity \f$\xi_w\f$
    double alpha;           ///< phase transition strength \f$\alpha_N\f$
    /// type of the hydrodynamic solution
    profile_calculator::HydroSolutionType sol_type; 

    /**
     * @brief Check consistency and validity of the fluid profile data.
     * @throws std::runtime_error if the check fails.
     */
    void check_profiles() const;

}; // class FluidProfile

} // namespace inverse_pt