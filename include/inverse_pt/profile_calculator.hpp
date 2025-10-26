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

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "inverse_pt/constants.hpp"
#include "inverse_pt/settings.hpp"

/**
 * @namespace profile_calculator
 * @brief Functions to calculate the fluid profile of a bubble in a cosmological
 *     phase transition.
 * 
 * This is based on 
 * 1. @anchor ref_Espinosa_2010
 *    J. R. Espinosa and J. M. No,
 *    "Energy Budget of Cosmological First-order Phase Transitions,"
 *    [JCAP **06** (2010), 028](https://doi.org/10.1088/1475-7516/2010/06/028),
 *    [arXiv:1004.4187 [hep-ph]](https://arxiv.org/abs/1004.4187).
 * 2. @anchor ref_Barni_2024
 *    G. Barni, S. Blasi and M. Vanvlasselaer, 
 *    "The hydrodynamics of inverse phase transitions,"
 *    [JCAP **10** (2024), 042](https://doi.org/10.1088/1475-7516/2024/10/042),
 *    [arXiv:2406.01596 [hep-ph]](https://arxiv.org/abs/2406.01596).
 * 3. @anchor ref_Barni_2025
 *    G. Barni, S. Blasi and M. Vanvlasselaer, 
 *    "Inverse bubbles from broken symmetries,"
 *    [arXiv:2503.01951 [hep-ph]](https://arxiv.org/abs/2503.01951).
 */
namespace inverse_pt::profile_calculator {

/// @brief Types of the hydrodynamic solutions.
enum class HydroSolutionType { 
    Deflagration,        ///< weak (subsonic) deflagration 
    Detonation,          ///< weak (supersonic) detonation 
    Hybrid,              ///< hybrid solution (supersonic deflagration)
    InverseDeflagration, ///< weak (supersonic) inverse deflagration
    InverseDetonation,   ///< weak (subsonic) inverse detonation
    InverseHybrid   ///< inverse hybrid solution (subsonic inverse deflagration)
};

/** 
 * @brief Selection which differential equation to solve when calculating the
 *     velocity profile.
 */
enum class ProfileSolveMode {
    D_XI_D_V,           ///< always solve \f$\frac{d\xi}{dv}\f$
    D_V_D_XI,           ///< always solve \f$\frac{dv}{d\xi}\f$
    /** 
     * solve \f$\frac{dv}{d\xi}\f$ when the solution is expected to end in a
     * shock and \f$\frac{d\xi}{dv}\f$ otherwise 
     */
    D_V_D_XI_FOR_SHOCKS 
};

/**
 * @brief Convert \c HydroSolutionType \c enum value to a human-readable string.
 * @param type The \c HydroSolutionType value
 * @return A <tt>const char*</tt> string representation of the enum value
 */
const std::string to_string(HydroSolutionType type);

/**
 * @brief Stream output operator for \c HydroSolutionType.
 * @param os The output stream
 * @param type The \c HydroSolutionType value to print
 * @return The output stream with the string representation of the
 *     \c HydroSolutionType appended
 */
std::ostream& operator<<(std::ostream& os, HydroSolutionType type);

/**
 * @brief Check if two numbers are equal.
 * 
 * This functions checks whether the two numbers \p x and \p y are equal within
 * the relative precision \p eps_rel and the absolute precision \p eps_abs.
 * 
 * @param x First number
 * @param y Second number
 * @param eps_rel Relative precision in the comparison (default: 1e-9)
 * @param eps_abs Absolute precision in the comparison (default: 1e-12)
 * @return Whether the two numbers are equal within the tolerance
 */
bool are_equal(double x, double y, double eps_rel=1e-9, double eps_abs=1e-12);

/**
 * @brief Compute the Lorentz boost factor \f$\gamma=\frac{1}{\sqrt{1-v^2}}\f$.
 * @param v Velocity \f$v\f$
 * @return Boost factor \f$\gamma\f$
 */
double gamma(double v);

/**
 * @brief Lorentz boosted velocity.
 * 
 * This boosts the velocity \f$v\f$ to the frame moving with velocity \f$\xi\f$,
 * Eq. (2.28) of @ref ref_Espinosa_2010 "Espinosa and No (2010)".
 * \f[ \mu(\xi,v) = \frac{\xi-v}{1-\xi v} \f]
 * 
 * @param xi Relative velocity \f$\xi\f$ of the new frame to which the velocity
 *           is boosted
 * @param v Velocity \f$v\f$ in the original frame
 * @return Boosted velocity \f$\mu(\xi, v)\f$
 */
double mu(double xi, double v);

/**
 * @brief Check if the velocity is between 0 and 1.
 * @param v Velocity to check
 * @return @c true if the velocity lies between 0 and 1, @c false otherwise
 */
bool check_velocity(double v);

/**
 * @brief Calculate the distance to the shock front.
 * 
 * This function returns the shock condition
 * \f$\mu(\xi_\text{sh}, v_\text{sh}) \xi_\text{sh} - c_s^2\f$.
 * It returns zero if \f$(\xi,v)\f$ corresponds to a shock front, a negative
 * value if the point is ahead of a shock front and a positive value if it is
 * behind the shock front.
 * 
 * @param xi Self-similar coordinate \f$\xi=r/t\f$
 * @param v Fluid vlocity \f$v\f$
 * @return Distance to the shock front
 */
double shock_condition(double xi, double v);

/**
 * @brief Calculate the distance to the Jouguet front.
 * 
 * This function returns the Jouguet condition
 * \f$\mu(\xi_\text{J}, v_\text{J})^2 - c_s^2\f$.
 * It returns zero if \f$(\xi,v)\f$ corresponds to a Jouguet front, a negative
 * value if the point is ahead of a Jouguet front and a positive value if it is
 * behind the Jouguet front.
 * 
 * @param xi Self-similar coordinate \f$\xi=r/t\f$
 * @param v Fluid vlocity \f$v\f$
 * @return Distance to the Jouguet front
 */
double jouguet_condition(double xi, double v);

/** 
 * @brief Calculate \f$v_+\f$ from \f$v_-\f$ and \f$\alpha_+\f$.
 * 
 * This function calculates the fluid velocity \f$v_+\f$ ahead of the bubble (or
 * shock) front from the velocity \f$v_-\f$ behind the front and the transition
 * strength \f$\alpha_+\f$ ahead of the front. This is Eq. (2.21) of 
 * @ref ref_Espinosa_2010 "Espinosa and No (2010)".
 * 
 * @param vm Fluid velocity \f$v_-\f$ behind the front
 * @param alpha Transition strength \f$\alpha_+\f$ ahead of the front
 * @return The solutions on the two branches for the fluid velocity \f$v_+\f$
 *     ahead of the front
 */
std::pair<double, double> get_v_plus_from_v_minus(double vm, double alpha);

/** 
 * @brief Calculate \f$v_-\f$ from \f$v_+\f$ and \f$\alpha_+\f$.
 * 
 * This function calculates the fluid velocity \f$v_-\f$ behind the bubble (or
 * shock) front from the velocity \f$v_+\f$ and the transition strength
 * \f$\alpha_+\f$ ahead of the front. This is Eq. (2.21) of
 * @ref ref_Espinosa_2010 "Espinosa and No (2010)" solved for \f$v_-\f$.
 * 
 * @param vp Fluid velocity \f$v_+\f$ ahead of the front
 * @param alpha Transition strength \f$\alpha_+\f$ ahead of the front
 * @return The solutions on the two branches for the fluid velocity \f$v_-\f$
 *     behind the front
 */
std::pair<double, double> get_v_minus_from_v_plus(double vp, double alpha);

/**
 * @brief Compute the transition strength from the fluid velocities in the wall
 *     rest frame.
 * 
 * Cf. Eq. (33) of @ref ref_Barni_2024 "Barni et al. (2024)".
 * 
 * @param vp Fluid velocity \f$v_+\f$ ahead of the wall
 * @param vm Fluid velocity \f$v_-\f$ behind the wall
 * @return Phase transition strength \f$\alpha_+\f$ ahead of the bubble
 */
double compute_alpha_plus(double vp, double vm);

/**
 * @brief Compute the Jouguet velocity for a direct phase transition, Eq. (36)
 *     of @ref ref_Barni_2024 "Barni et al. (2024)".
 * @param alpha Strength parameter \f$\alpha_+\f$
 * @return Jouguet velocity \f$v_J^\text{direct}\f$
 */
double direct_jouguet_velocity(double alpha);

/**
 * @brief Compute the Jouguet velocity for an inverse phase transition, solution
 *     of Eq. (39) of @ref ref_Barni_2024 "Barni et al. (2024)".
 * @param alpha Strength parameter \f$\alpha_+\f$
 * @return Inverse Jouguet velocity \f$v_J^\text{inverse}\f$
 */
double inverse_jouguet_velocity(double alpha);

/** 
 * @brief Compute the maximal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     detonation solutions.
 * 
 * This functions returns the maximal value of the phase transitions strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for direct detonations.
 * 
 * For direct detonations \f$\alpha_+\f$ and \f$\alpha_N\f$ are identical, and
 * the transition strength must be smaller than the (direct) Jouguet strength,
 * \f$\alpha_{N/+} < \alpha_J^\text{direct}(\xi_w) = \alpha_+(\xi_w, c_s)\f$.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_N^\mathrm{max}\f$ (\c true, default) or
 *     \f$\alpha_+^\mathrm{max}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-4)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return The maximal value of the transition strength for direct detonations
 */
double alpha_max_detonation(
    double xi_w, 
    bool return_alpha_N=true, 
    double step_size=1e-4,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/** 
 * @brief Compute the maximal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     deflagration solutions.
 * 
 * This functions returns the maximal value of the phase transitions strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for direct deflagrations.
 * 
 * For \f$\alpha_+ > \frac{1}{3}\f$, the subsonic branch of the fluid velocity
 * \f$v_+(\xi_w, \alpha_+)\f$ ahead of the wall becomes negative.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_N^\mathrm{max}\f$ (\c true, default) or
 *     \f$\alpha_+^\mathrm{max}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-4)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return The maximal value of the transition strength for direct deflagrations
 */
double alpha_max_deflagration(
    double xi_w, 
    bool return_alpha_N=true, 
    double step_size=1e-4,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/** 
 * @brief Compute the maximal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     hybrid solutions.
 * 
 * This functions returns the maximal value of the phase transitions strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is 
 * \c false) for direct hybrids.
 * 
 * For \f$\alpha_+ > \frac{1}{3}\f$, the subsonic branch of the fluid velocity
 * \f$v_+(c_s, \alpha_+)\f$ ahead of the wall becomes negative.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_N^\mathrm{max}\f$ (\c true, default) or
 *     \f$\alpha_+^\mathrm{max}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-4)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return The maximal value of the transition strength for direct hybrids
 */
double alpha_max_hybrid(
    double xi_w, 
    bool return_alpha_N=true, 
    double step_size=1e-4,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/** 
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     hybrid solutions.
 * 
 * This functions returns the minimal value of the phase transitions strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for direct hybrids.
 * 
 * In order to be able to evolve towards a shock front, the matching solution at
 * the bubble wall needs to satisfy \f$\xi_w v_+(c_s, \alpha_+) < c_s^2\f$,
 * hence \f$\alpha_+ > \frac{\left(\xi_w - c_s\right)^2}{3 \xi_w^2 - c_s^2}\f$.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_N^\mathrm{min}\f$ (\c true, default) or
 *     \f$\alpha_+^\mathrm{min}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-4)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return The minimal value of the transition strength for direct hybrids
 */
double alpha_min_hybrid(
    double xi_w, 
    bool return_alpha_N=true,
    double step_size=1e-4,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/** 
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     inverse deflagration solutions.
 * 
 * This functions returns the minimal value of the phase transitions strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for inverse deflagrations.
 * 
 * For inverse deflagrations \f$\alpha_+\f$ and \f$\alpha_N\f$ are identical.
 * For \f$\alpha_+ < \frac{c_s^2 - xi_w}{1 + \xi_w}\f$, the supersonic branch of
 * the fluid velocity \f$v_-(\xi_w, \alpha_+)\f$ behind the wall becomes
 * superluminal.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_N^\mathrm{min}\f$ (\c true, default) or
 *     \f$\alpha_+^\mathrm{min}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-4)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return The minimal value of the transition strength for inverse 
 *     deflagrations
 */
double alpha_min_inverse_deflagration(
    double xi_w, 
    bool return_alpha_N=true, 
    double step_size=1e-4,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/** 
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     inverse hybrid solutions.
 * 
 * This functions returns the minimal value of the phase transitions strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for inverse hybrids.
 * 
 * For \f$\alpha_+ < 1 - 2 c_s\f$, the supersonic branch of the fluid velocity
 * \f$v_-(c_s, \alpha_+)\f$ behind the wall becomes superluminal.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_N^\mathrm{min}\f$ (\c true, default) or
 *     \f$\alpha_+^\mathrm{min}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-4)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return The minimal value of the transition strength for inverse hybrids
 */
double alpha_min_inverse_hybrid(
    double xi_w, 
    bool return_alpha_N=true, 
    double step_size=1e-4,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/** 
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for
 *     inverse detonation solutions.
 * 
 * This functions returns the minimal value of the phase transitions strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is 
 * \c false) for inverse detonations.
 * 
 * For inverse detonations, the transition strength must be larger than the
 * inverse Jouguet strength,
 * \f$\alpha_+ > \alpha_J^\text{inv} = \alpha_+(c_s, \xi_w)\f$.
 * 
 * @note This is also the maximal strength for inverse hybrids.
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param return_alpha_N Whether the function should return 
 *     \f$\alpha_N^\mathrm{min}\f$ (\c true, default) or
 *     \f$\alpha_+^\mathrm{min}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-4)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return The minimal value of the transition strength for inverse detonations
 */
double alpha_min_inverse_detonation(
    double xi_w, 
    bool return_alpha_N=true, 
    double step_size=1e-4,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/** 
 * @brief Compute the maximal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for the
 *     given solution type.
 * 
 * This functions returns the maximal value of the phase transitions strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for the solution type specified in @p sol_type.
 * 
 * @param sol_type Type of the solution
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_N^\mathrm{max}\f$ (\c true, default) or
 *     \f$\alpha_+^\mathrm{max}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-4)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return The maximal value of the transition strength for the specified
 *     solution type
 */
double alpha_max(
    HydroSolutionType sol_type, 
    double xi_w, 
    bool return_alpha_N=true, 
    double step_size=1e-4,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/** 
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$ for the
 *     given solution type.
 * 
 * This functions returns the minimal value of the phase transitions strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false) for the solution type specified in @p sol_type.
 *
 * @param sol_type Type of the solution
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_N^\mathrm{min}\f$ (\c true, default) or
 *     \f$\alpha_+^\mathrm{min}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-4)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return The minimal value of the transition strength for the specified
 *     solution type
 */
double alpha_min(
    HydroSolutionType sol_type,
    double xi_w,
    bool return_alpha_N=true,
    double step_size=1e-4,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/** 
 * @brief Compute the maximal value of \f$\alpha_N\f$ or \f$\alpha_+\f$.
 * 
 * This functions returns the maximal value of the phase transitions strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is
 * \c false).
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_N^\mathrm{max}\f$ (\c true, default) or
 *     \f$\alpha_+^\mathrm{max}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-4)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return The maximal value of the transition strength
 */
double alpha_max(
    double xi_w, 
    bool return_alpha_N=true, 
    double step_size=1e-4,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/** 
 * @brief Compute the minimal value of \f$\alpha_N\f$ or \f$\alpha_+\f$.
 * 
 * This functions returns the minimal value of the phase transitions strength
 * \f$\alpha_N\f$ far ahead of the wall (if @p return_alpha_N is \c true) or
 * \f$\alpha_+\f$ directly in front of the wall (if @p return_alpha_N is 
 * \c false).
 * 
 * @param xi_w Wall velocity \f$\xi_w\f$
 * @param return_alpha_N Whether the function should return
 *     \f$\alpha_N^\mathrm{min}\f$ (\c true, default) or
 *     \f$\alpha_+^\mathrm{min}\f$ (\c false)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-4)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return The minimal value of the transition strength
 */
double alpha_min(
    double xi_w, 
    bool return_alpha_N=true, 
    double step_size=1e-4,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/**
 * @brief Classify the type of the solution based on \f$\xi_w\f$ and
 *     \f$\alpha_N\f$ in the fluid frame.
 * @param xi_wall Bubble wall velocity \f$\xi_w\f$
 * @param alpha Phase transition strength \f$\alpha_N\f$ (far ahead of the wall)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-4)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default: 
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return Type of the hydrodynamic solution 
 */
HydroSolutionType classify_wave_fluid_frame(
    double xi_wall, 
    double alpha, 
    double step_size=1e-4,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/**
 * @brief Classify the type of the solution based on \f$v_+\f$ and \f$v_-\f$ in
 *     the wall rest frame.
 * @param vp Fluid velocity \f$v_+\f$ ahead of the wall
 * @param vm Fluid velocity \f$v_-\f$ behind the wall
 * @return Type of the hydrodynamic solution 
 */
HydroSolutionType classify_wave_wall_frame(double vp, double vm);

/**
 * @brief Compute the fluid velocities \f$v_+\f$ and \f$v_-\f$ in the bubble
 *     rest frame.
 * 
 * @note For direct and inverse hybrids, \f$v_-\f$ (direct) and \f$v_+\f$
 *     (inverse) are fixed to $c_s$, respectively. This function then returns
 *     \f$\xi_w\f$ instead of \f$c_s\f$.
 * 
 * @param xi_wall Bubble wall velocity \f$\xi_w\f$
 * @param alpha Phase transition strength \f$\alpha_+\f$ ahead of the wall
 * @param sol_type Type of the solution.
 * @return The fluid velocities ahead of and behind the wall in the fluid rest
 *     frame, \f$v_+\f$ and \f$v_-\f$
 */
std::pair<double, double> compute_fluid_velocities(
    double xi_wall, double alpha, HydroSolutionType sol_type
);

/**
 * @brief Derivative \f$\frac{d\xi}{dv}\f$.
 * 
 * We solve this differential equation to obtain \f$\xi(v)\f$ and then invert to
 * obtain \f$v(\xi)\f$.
 * This is Eq. (2.27) of @ref ref_Espinosa_2010 "Espinosa and No (2010)".
 * 
 * @param v Fluid velocity \f$v\f$
 * @param xi self-similar coordinate \f$\xi=r/t\f$ 
 * @return Derivative \f$\frac{d\xi}{dv}\f$
 */
double dxi_dv(double v, double xi);

/**
 * @brief Derivative \f$\frac{dv}{d\xi}\f$.
 * 
 * We solve this differential equation to directly obtain \f$v(\xi)\f$.
 * This is Eq. (2.27) of @ref ref_Espinosa_2010 "Espinosa and No (2010)".
 * 
 * @param xi self-similar coordinate \f$\xi=r/t\f$ 
 * @param v Fluid velocity \f$v\f$
 * @return Derivative \f$\frac{dv}{d\xi}\f$
 */
double dv_dxi(double xi, double v);

/**
 * @brief Solve the differential equation \f$\frac{d\xi}{d v}\f$ for the
 *     velocity profile.
 * 
 * This function solves the differential equation using a fourth-order
 * Runge-Kutta method. The first two arguments are references to the vectors in
 * which the solutions are stored. These vectors must already contain the
 * initial condition.
 * 
 * @param xi_vec Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v_vec Reference to the vector containing the fluid velocity \f$v\f$
 * @param v_stop Final fluid velocity \f$v_f\f$
 * @param dv Fluid velocity step size
 * @param stop_at_shock Whether to stop the evolution when a shock front is
 *     encoutered (default: \c true)
 * @return A status code indicating the reason for termination:
 *     0  - reached v_stop
 *     1  - encountered shock 
 *     -1 - entered unphysical range for \f$\xi\f$
 *     -2 - reached super-luminal \f$v\f$
 */
int integrate_xi_of_v(
    std::vector<double>& xi_vec, 
    std::vector<double>& v_vec, 
    double v_stop, 
    double dv, 
    bool stop_at_shock=true
);

/**
 * @brief Solve the differential equation \f$\frac{dv}{d \xi}\f$ for the
 *     velocity profile.
 * 
 * This function solves the differential equation using a fourth-order
 * Runge-Kutta method. The first two arguments are references to the vectors in
 * which the solutions are stored. These vectors must already contain the
 * initial condition.
 * 
 * @param xi_vec Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v_vec Reference to the vector containing the fluid velocity \f$v\f$
 * @param v_stop Final fluid velocity \f$v_f\f$
 * @param dxi Self-similar coordinate step size
 * @param stop_at_shock Whether to stop the evolution when a shock front is
 *     encoutered (default: \c true)
 * @return A status code indicating the reason for termination:
 *     0  - reached v_stop
 *     1  - encountered shock 
 *     -1 - entered unphysical range for \f$\xi\f$
 *     -2 - reached super-luminal fluid velocity \f$v\f$
 *     -3 - crossed Jouguet front
 */
int integrate_v_of_xi(
    std::vector<double>& xi_vec,
    std::vector<double>& v_vec,
    double v_stop,
    double dxi,
    bool stop_at_shock=true
);

/**
 * @brief Solve the differential equation for the velocity profile.
 * 
 * This function solves the differential equation using a fourth-order
 * Runge-Kutta method. The first two arguments are references to the vectors in
 * which the solutions are stored. These vectors must already contain the
 * initial condition.
 * 
 * @param xi_vec Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v_vec Reference to the vector containing the fluid velocity \f$v\f$
 * @param v_stop Final fluid velocity \f$v_f\f$
 * @param step_size Fluid velocity or self-similar coordinate step size
 *     (depending on @p solve_dxi_dv)
 * @param stop_at_shock Whether to stop the evolution when a shock front is
 *     encoutered (default: \c true)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c ProfileSolveMode::D_XI_D_V)
 * @return A status code indicating the reason for termination:
 *     0  - reached v_stop
 *     1  - encountered shock 
 *     -1 - entered unphysical range for \f$\xi\f$
 *     -2 - reached super-luminal fluid velocity \f$v\f$
 *     -3 - crossed Jouguet front
 */
int integrate_velocity_profile(
    std::vector<double>& xi_vec,
    std::vector<double>& v_vec,
    double v_stop,
    double step_size,
    bool stop_at_shock=true,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_XI_D_V
);

/**
 * @brief Compute the velocity profile for a direct detonation.
 * 
 * This functions solves the fluid velocity profile for a direct detonation and
 * stores the results of the velocity \f$v\f$ as a function of the self-similar
 * coordiate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param v_plus Fluid velocity \f$v_+\f$ in front of the wall
 * @param v_minus Fluid velocity \f$v_-\f$ behind the wall
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-3)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return Index of the first element outside the bubble (i.e. at the wall
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile_direct_detonation(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double v_minus, 
    double step_size=1e-3,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/**
 * @brief Compute the velocity profile for a direct deflagration.
 * 
 * This functions solves the fluid velocity profile for a direct deflagration
 * and stores the results of the velocity \f$v\f$ as a function of the 
 * self-similar coordiate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param v_plus Fluid velocity \f$v_+\f$ in front of the wall
 * @param v_minus Fluid velocity \f$v_-\f$ behind the wall
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-3)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return Index of the first element outside the bubble (i.e. at the wall
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile_direct_deflagration(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double v_minus, 
    double step_size=1e-3,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/**
 * @brief Compute the velocity profile for a direct hybrid.
 * 
 * This functions solves the fluid velocity profile for a direct hybrid and
 * stores the results of the velocity \f$v\f$ as a function of the self-similar
 * coordiate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param v_plus Fluid velocity \f$v_+\f$ in front of the wall
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-3)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return Index of the first element outside the bubble (i.e. at the wall
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile_direct_hybrid(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double xi_wall, 
    double step_size=1e-3,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/**
 * @brief Compute the velocity profile for an inverse detonation.
 * 
 * This functions solves the fluid velocity profile for an inverse detonation
 * and stores the results of the velocity \f$v\f$ as a function of the
 * self-similar coordiate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param v_plus Fluid velocity \f$v_+\f$ in front of the wall
 * @param v_minus  Fluid velocity \f$v_-\f$ behind the wall
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-3)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return Index of the first element outside the bubble (i.e. at the wall
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile_inverse_detonation(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double v_minus, 
    double step_size=1e-3,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/**
 * @brief Compute the velocity profile for an inverse deflagration.
 * 
 * This functions solves the fluid velocity profile for an inverse deflagration
 * and stores the results of the velocity \f$v\f$ as a function of the
 * self-similar coordiate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param v_plus Fluid velocity \f$v_+\f$ in front of the wall
 * @param v_minus Fluid velocity \f$v_-\f$ behind the wall
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-3)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return Index of the first element outside the bubble (i.e. at the wall 
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile_inverse_deflagration(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double v_plus, 
    double v_minus, 
    double step_size=1e-3,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/**
 * @brief Compute the velocity profile for an inverse hybrid.
 * 
 * This functions solves the fluid velocity profile for an inverse hybrid and
 * stores the results of the velocity \f$v\f$ as a function of the self-similar
 * coordiate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param v_minus Fluid velocity \f$v_-\f$ behind the wall
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-3)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return Index of the first element outside the bubble (i.e. at the wall
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile_inverse_hybrid(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double xi_wall, 
    double v_minus, 
    double step_size=1e-3,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/**
 * @brief Get the position of the shock front.
 * 
 * @param sol_type Type of the solution.
 * @param size Length of the profile vectors
 * @return Index of the first element outside the shock front
 */
size_t get_shock_index(HydroSolutionType sol_type, size_t size);

/**
 * @brief Get the coordinate and velocity of the shock front.
 * 
 * @param sol_type Type of the solution.
 * @param xi Self-similar coordinate \f$\xi=r/t\f$ of the profile
 * @param v Fluid velocity profile \f$v(\xi)\f$
 * @return The shock front position \f$(\xi_\text{sh}, v_\text{sh})\f$
 */
std::pair<double, double> get_shock_front(
    HydroSolutionType sol_type, 
    const std::vector<double>& xi, 
    const std::vector<double>& v
);

/**
 * @brief Compute the velocity profile for a bubble in a cosmological phase
 *     transition.
 * 
 * This functions solves the fluid velocity profile for a cosmological phase
 * transition and stores the results of the velocity \f$v\f$ as a function of
 * the self-similar coordiate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param alpha_N Transition strength \f$\alpha_N\f$ (far in front of the wall)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-3)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return Index of the first element outside the bubble (i.e. at the wall 
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile(
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double xi_wall, 
    double alpha_N, 
    double step_size=1e-3,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/**
 * @brief Compute the velocity profile for a bubble in a cosmological phase
 *     transition.
 * 
 * This functions solves the fluid velocity profile for a cosmological phase
 * transition and stores the results of the velocity \f$v\f$ as a function of
 * the self-similar coordiate \f$\xi\f$ in the vectors @p xi and @p v.
 * 
 * @param sol_type Type of the solution.
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param alpha_p Transition strength \f$\alpha_+\f$ ahead of the wall
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-3)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return Index of the first element outside the bubble (i.e. at the wall 
 *     velocity \f$\xi_w\f$)
 */
size_t compute_velocity_profile(
    HydroSolutionType sol_type, 
    std::vector<double>& xi, 
    std::vector<double>& v, 
    double xi_wall, 
    double alpha_p, 
    double step_size=1e-3,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/** 
 * @brief Calculate a segment of the enthalpy density integral.
 * 
 * This function calculates the segment of the integral in Eq. (2.29) of
 * @ref ref_Espinosa_2010 "Espinosa and No (2010)"
 * \f[
 *  I = \left(1-\frac{1}{c_s^2}\right) \int\limits_{v_1}^{v_2} \gamma^2 \mu dv 
 *    = \left(1-\frac{1}{c_s^2}\right) \int\limits_{v_1}^{v_2} \frac{1}{1-v^2}
 *      \frac{\xi - v}{1 - \xi v} \mu dv
 * \f]
 * assuming that \f$v\f$ is linear in the segment \f$(\xi_1,\xi_2)\f$, i.e.
 * \f$\xi = \xi_0 + \partial_v \xi\, v\f$, where
 * \f$\xi_0=\frac{\xi_1 v_2 - \xi_2 v_1}{v_2-v_1}\f$ and
 * \f$\partial_v \xi = \frac{\xi_2-\xi_1}{v_2-v_1}\f$.
 * The integral can the be calculated using partial fraction decomposition:
 * \f[
 *  I = \left(1-\frac{1}{c_s^2}\right) \int\limits_{v_1}^{v_2} \left[ 
 *     \frac{v}{1-v^2} 
 *     + \frac{\partial_v \xi\,v + \xi_0}{\partial_v \xi\, v^2 + \xi_0 v - 1} 
 *   \right] dv
 *    =  \left(1-\frac{1}{c_s^2}\right) \left\{ 
 *          \frac{1}{2} \int\limits_{1-v_2^2}^{1-v_1^2} \frac{du}{u} 
 *          + \frac{1}{2} \int\limits_{1-\xi_2 v_2}^{1-\xi_1 v_1} \frac{du}{u}
 *          + \int\limits_{v_1}^{v_2} 
 *              \frac{\xi_0}{\partial_v \xi\, v^2 + \xi_0 v - 1} dv
 *       \right\} .
 * \f]
 * 
 * @param xi1 Self-similar coordinate \f$\xi_1\f$ at the beginning of the 
 *     segment
 * @param xi2 Self-similar coordinate \f$\xi_2\f$ at the end of the segment
 * @param v1 Fluid velocity \f$v_1\f$ at the beginning of the segment
 * @param v2 Fluid velocity coordinate \f$v_2\f$ at the end of the segment
 */
double integrate_enthalpy_segment(double xi1, double xi2, double v1, double v2);

/**
 * @brief Compute the enthalpy density profile.
 * 
 * Compute the enthalpy density profile from the fluid velocity profile,
 * Eq. (2.29) of @ref ref_Espinosa_2010 "Espinosa and No (2010)".
 * 
 * @param xi Reference to the vector containing the self-similar coordinate 
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$
 * @param idx_wall Index of the first element outside the bubble
 * @param idx_shock Index of the first element outside the shock front (use 
 *     \c constants::NO_INDEX if there is no shock front)
 * @return Enthalpy density profile \f$w(\xi)\f$
 */
std::vector<double> compute_enthalpy_density(
    const std::vector<double>& xi, 
    const std::vector<double>& v, 
    size_t idx_wall, 
    size_t idx_shock=constants::NO_INDEX
);

/**
 * @brief Compute the energy density profile.
 * 
 * We assume the Bag equation of state, i.e. \f$e = \frac{3}{4} w + \epsilon\f$.
 * 
 * @param w Reference to the vector containing the enthalpy density \f$w\f$
 * @param idx_wall Index of the first element outside the bubble
 * @param epsilon Bag constant
 * @return Energy density profile \f$e(\xi)\f$
 */
std::vector<double> compute_energy_density_profile(
    const std::vector<double>& w, size_t idx_wall, double epsilon
);

/**
 * @brief Compute the pressure density profile.
 * 
 * We assume the Bag equation of state, i.e. \f$p = \frac{1}{4} w - \epsilon\f$.
 * 
 * @param w Reference to the vector containing the enthalpy density \f$w\f$
 * @param idx_wall Index of the first element outside the bubble
 * @param epsilon Bag constant
 * @return Pressure density profile \f$p(\xi)\f$
 */
std::vector<double> compute_pressure_density_profile(
    const std::vector<double>& w, size_t idx_wall, double epsilon
);

/**
 * @brief Compute the Bag constant (vacuum energy density) 
 *     \f$\epsilon = \frac{3}{4} \alpha_+ w_+ = \frac{3}{4} \alpha_N w_N\f$.
 * 
 * @param w Enthalpy density, either \f$w_+\f$ (directly in front) or \f$w_N\f$
 *     (far in front of the bubble wall)
 * @param alpha Transition strength, either \f$\alpha_+\f$ (directly in front)
 *     or \f$\alpha_N\f$ (far in front of the bubble wall)
 * @return Bag constant \f$\epsilon\f$
 */
double compute_epsilon(double w, double alpha);

/**
 * @brief Compute the transition strength \f$\alpha_N\f$ from the value in front
 *     of the bubble wall \f$\alpha_+\f$.
 * 
 * @param sol_type Type of the solution.
 * @param alpha_plus Transition strength \f$\alpha_+\f$ immediately in front of
 *     the bubble wall
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return Transition strength \f$\alpha_N\f$ far in front of the bubble wall
 */
double compute_alpha_N(
    HydroSolutionType sol_type, 
    double alpha_plus, 
    double xi_wall, 
    double step_size,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/**
 * @brief Find the transition strength \f$\alpha_+\f$ from the strength
 *     \f$\alpha_N\f$ in the false vacuum phase.
 * 
 * If the solution type is not provided, the function calls 
 * @ref classify_wave_fluid_frame internally.
 * 
 * @param alpha_N Transition strength \f$\alpha_N\f$ far in front of the wall
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 */
double find_alpha_plus(
    double alpha_N, 
    double xi_wall, 
    double step_size,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/**
 * @brief Find the transition strength \f$\alpha_+\f$ from the strength
 *     \f$\alpha_N\f$ in the false vacuum phase.
 * 
 * @param sol_type Type of the solution.
 * @param alpha_N Transition strength \f$\alpha_N\f$ far in front of the wall
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return Value of \f$\alpha_+\f$
 */
double find_alpha_plus(
    HydroSolutionType sol_type,
    double alpha_N, 
    double xi_wall, 
    double step_size,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/**
 * @brief Compute the efficiency factor \f$\kappa_\mathrm{sw}\f$.
 * 
 * Compute the efficiency factor \f$\kappa_\mathrm{sw}\f$ for converting vacuum
 * energy to bulk kinetic energy, Eq. (30) of
 * @ref ref_Espinosa_2010 "Espinosa and No (2010)" or Eq. (49) of
 * @ref ref_Barni_2024 "Barni et al. (2024)".
 * 
 * @param xi Reference to the vector containing the self-similar coordinate
 *     \f$\xi\f$
 * @param v Reference to the vector containing the fluid velocity \f$v\f$ 
 * @param w Reference to the vector containing the enthalpy density \f$w\f$
 * @param idx_wall Index of the first element outside the bubble
 * @param epsilon Bag constant
 * @return Efficiency factor \f$\kappa_\mathrm{sw}\f$
 */
double compute_efficiency_factor(
    const std::vector<double>& xi, 
    const std::vector<double>& v, 
    const std::vector<double>& w, 
    size_t idx_wall,
    double epsilon
);

/**
 * @brief Compute the ratio \f$\Psi\f$ of degrees of freedom in local thermal
 *     equilibrium.
 * 
 * Compute the ration \f$\Psi = \frac{a_-}{a_+}\f$ of degrees of freedom ahead
 * of and behind the wall. This assumes local thermal equilibrium and uses
 * entropy conservation to relate the temperatures on both sides of the wall,
 * i.e. \f$s_+ v_+ \gamma_+ = s_- v_- \gamma_-\f$.
 * 
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param alpha_N Transition strength \f$\alpha_N\f$ (far in front of the wall)
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-4)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return Ratio \f$\Psi\f$ of degrees of freedom across the wall
 */
double compute_dofs_ratio_LTE(
    double xi_wall, 
    double alpha_N, 
    double step_size=1e-4,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/**
 * @brief Find the LTE wall velocity \f$\xi_w\f$ for a given dof ratio
 *     \f$\Psi\f$.
 * 
 * Find the wall velocity \f$\xi_w\f$ for a given ratio of degrees of freedom
 * \f$\Psi = \frac{a_-}{a_+}\f$ assuming local thermal equilibrium.
 * 
 * @param alpha_N Transition strength \f$\alpha_N\f$ (far in front of the wall)
 * @param Psi Ratio of degrees of freedom \f$\Psi\f$
 * @param xi_min Minimal wall velocity \f$\xi_w^\mathrm{min}\f$ to consider
 * @param xi_max Maximal wall velocity \f$\xi_w^\mathrm{max}\f$ to consider
 * @param step_size Step size for the fluid velocity (or coordinate) when
 *     solving the differential equation (default: 1e-4)
 * @param solve_mode \c ProfileSolveMode switch determining which differential
 *     equation for the velocity profile is solved (default:
 *     \c ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
 * @return Wall velocity \f$\xi_w\f$
 */
double find_LTE_wall_velocity(
    double alpha_N, 
    double Psi, 
    double xi_min, 
    double xi_max, 
    double step_size=1e-4,
    ProfileSolveMode solve_mode=ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
);

/** 
 * @brief Check profile for consistency.
 * 
 * This function returns @c true if @p xi and @p profile have the same size and
 * at least size 6. Otherwise it throws an error (if @p throw_error is @c true)
 * or returns @c false.
 * 
 * @param xi Self-similar coordinate \f$\xi=r/t\f$ of the profile
 * @param profile Vector containing the fluid profile
 * @param throw_error Whether to throw an error if the test fails 
 * @return Whether the profile passes the test
 */
bool check_profiles(
    const std::vector<double>& xi, 
    const std::vector<double>& profile, 
    bool throw_error=true
);


} // namespace profile_calculator