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

#include "inverse_pt/fluid_profile.hpp"
#include "inverse_pt/sound_shell_model.hpp"

#include <utility>
#include <vector>

namespace inverse_pt {

/**
 * @brief Class to calculate the velocity and graviational wave spectrum in the
 *     sound shell model.
 */
class SoundShellSpectrum {
public:

    /// @name Constructors
    ///@{
    
    /**
     * @brief Constructor for the @ref SoundShellSpectrum from a
     *     @ref FluidProfile instance.
     * @param profile @ref FluidProfile instance
     * @param R_H Mean bubble radius at collision times Hubble rate,
     *     \f$H_* R_*\f$
     * @param nhist Nucleation history (default:
     *     \c sound_shell_model::NucleationHistory::EXPONENTIAL)
     * @param exhist Expansion history (default:
     *     \c sound_shell_model::ExpansionHistory::RADIATION_DOMINATION)
     */
    SoundShellSpectrum(
        const FluidProfile& profile, 
        double R_H,
        sound_shell_model::NucleationHistory nhist 
            = sound_shell_model::NucleationHistory::EXPONENTIAL,
        sound_shell_model::ExpansionHistory exhist 
            = sound_shell_model::ExpansionHistory::RADIATION_DOMINATION
    );

    /** 
     * @brief Constructor for the @ref SoundShellSpectrum from the wall velocity
     *     \f$\xi_w\f$ and strength parameter \f$\alpha_N\f$.
     * @param v_w wall veclocity \f$\xi_w\f$
     * @param alpha_N strength parameter \f$\alpha_N\f$ (far in front of the
     *     wall)
     * @param R_H Mean bubble radius at collision times Hubble rate,
     *     \f$H_* R_*\f$
     * @param nhist Nucleation history (default: 
     *     \c sound_shell_model::NucleationHistory::exponential)
     * @param exhist Expansion history (default: 
     *     \c sound_shell_model::ExpansionHistory::RADIATION_DOMINATION)
     * @param step_size Step size for the fluid velocity (or coordinate) when
     *     solving the differential equation (default: 1e-3)
     * @param solve_mode \c profile_calculator::ProfileSolveMode switch
     *     determining which differential equation for the velocity profile is
     *     solved (default: 
     *     \c profile_calculator::ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)     
     */
    SoundShellSpectrum(
        double v_w, 
        double alpha_N, 
        double R_H,
        sound_shell_model::NucleationHistory nhist 
            = sound_shell_model::NucleationHistory::EXPONENTIAL,
        sound_shell_model::ExpansionHistory exhist 
            = sound_shell_model::ExpansionHistory::RADIATION_DOMINATION,
        double step_size=1e-3,
        profile_calculator::ProfileSolveMode solve_mode
            = profile_calculator::ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
    );

    ///@}

    /**
     * @brief Set the values of \f$\chi = k T\f$ at which
     *     \f$|\mathcal{A}_+(\chi)|^2\f$ is evaluated.
     * 
     * The values of \f$\chi = k T\f$ are spaced logarithmically in @p N steps
     * between (the logarithms of) @p chi_min and @p chi_max.
     * 
     * @param chi_min Minimal value of \f$\chi = k T\f$
     * @param chi_max Maximal value of \f$\chi = k T\f$
     * @param N Number of points
     */
    void set_chi(double chi_min, double chi_max, size_t N);

    /**
     * @brief Set the values of \f$\chi = k T\f$ at which
     *     \f$|\mathcal{A}_+(\chi)|^2\f$ is evaluated.
     * 
     * @param chi_vec Vector with values of \f$\chi = k T\f$
     */
    void set_chi(const std::vector<double>& chi_vec);

    /**
     * @brief Set the values of the dimensionless momentum \f$q R_*\f$ at which
     *     the kinetic spectrum is evaluated.
     * 
     * The values of \f$q R_*\f$ are spaced logarithmically in @p N steps
     * between (the logarithms of) @p qR_min and @p qR_max.
     * 
     * @param qR_min Minimal value of \f$q R_*\f$
     * @param qR_max Maximal value of \f$q R_*\f$
     * @param N Number of points
     */
    void set_qR(double qR_min, double qR_max, size_t N);

    /**
     * @brief Set the values of the dimensionless momentum \f$q R_*\f$ at which
     *     the kinetic spectrum is evaluated.
     * 
     * @param qR_vec Vector with values of \f$q R_*\f$
     */
    void set_qR(const std::vector<double>&qR_vec);

    /**
     * @brief Set the values of the dimensionless momentum \f$q R_*\f$ at which
     *     the kinetic spectrum is evaluated based on the values of \f$k R_*\f$
     *     at which the graviational wave spectrum is calculated.
     * 
     * This function sets \f$q R_*\f$ making sure that the entire \f$q R_*\f$
     * range required to calculate the gravitational wave spectrum is covered.
     * It sets \f$q R_*\f$ in @p N logarithmically spaced points in the range
     * \f$(R_* q_\mathrm{low}, R_* q_\mathrm{high}\f$ with
     * \f[
     *   q_\substack{\mathrm{low}\\\mathrm{high}} = 
     *   k_\substack{\mathrm{low}\\\mathrm{high}}\,
     *   \frac{1 \mp c_s}{2\,c_s}\,,
     * \f]
     * where \f$R_* k_\mathrm{low}\f$ and \f$R_* k_\mathrm{high}\f$ lie one
     * logarithmic step blow and above the lowest and highest values of
     * \f$k R_*\f$, with the step size determined by the spacing between the 
     * first two and last two values, respectively.
     * 
     * @param N Number of points
     */
    void set_qR_from_kR(size_t N);

    /**
     * @brief Set the values of the dimensionless momentum \f$k R_*\f$ at which
     * the gravitational wave spectrum is evaluated.
     * 
     * The values of \f$k R_*\f$ are spaced logarithmically in @p N steps
     * between (the logarithms of) @p kR_min and @p kR_max. If @p N_qR is
     * provided, the values for \f$q R_*\f$ for the evaluation of the velocity
     * spectrum are set accordingly with @p N_qR points using
     * @ref set_qR_from_kR.
     * 
     * @param kR_min Minimal value of \f$k R_*\f$
     * @param kR_max Maximal value of \f$k R_*\f$
     * @param N Number of points (for \f$k R_*\f$)
     * @param N_qR If provided, @ref set_qR_from_kR is called to set \f$q R_*\f$
     *     with @p N_qR points
     */
    void set_kR(
        double kR_min, 
        double kR_max,
        size_t N, 
        size_t N_qR=static_cast<size_t>(-1)
    );

    /**
     * @brief Set the values of the dimensionless momentum \f$k R_*\f$ at which
     *     the gravitational wave spectrum is evaluated.
     * 
     * If @p N_qR is provided, the values for \f$q R_*\f$ for the evaluation of
     * the velocity spectrum are set accordingly with @p N_qR points using
     * @ref set_qR_from_kR.
     * 
     * @param kR_vec Vector with values of \f$k R_*\f$
     * @param N_qR If provided, @ref set_qR_from_kR is called to set \f$q R_*\f$
     *     with @p N_qR points
     */
    void set_kR(
        const std::vector<double>& kR_vec, size_t N_qR=static_cast<size_t>(-1)
    );

    /**
     * @brief Calculate the shape function \f$|\mathcal{A}_+(\chi)|^2\f$,
     *     Eq. (35) of 
     *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     * 
     * The shape function is evaluated on @p N logarithmically spaced values
     * between @p chi_min and @p chi_max. The values of \f$\chi = k T\f$ are
     * stored internally and can be retrieved calling @ref get_chi. The values
     * of \f$|\mathcal{A}_+(\chi)|^2\f$ are stored internally and can be
     * retrieved calling @ref get_shape_function. The function returns the
     * vectors containing the values of the functions \f$f'(\chi)\f$ and
     * \f$l(\chi)\f$ defined in Eqs. (30) and (31) of
     * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     * 
     * @param chi_min Minimal value of \f$\chi = k T\f$
     * @param chi_max Maximal value of \f$\chi = k T\f$
     * @param N Number of points
     * @return Pair of vectors containing the values of \f$f'(\chi)\f$ and
     *     \f$l(\chi)\f$
     */
    std::pair<std::vector<double>,std::vector<double>> calculate_A_squared(
        double chi_min, double chi_max, size_t N
    );

    /**
     * @brief Calculate the shape function \f$|\mathcal{A}_+(\chi)|^2\f$,
     *     Eq. (35) of
     *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     * 
     * The shape function is evaluated for the values provided in @p chi_vec.
     * The values of \f$\chi = k T\f$ are stored internally and can be retrieved
     * calling @ref get_chi. The values of \f$|\mathcal{A}_+(\chi)|^2\f$ are
     * stored interval and can be retrieved calling @ref get_shape_function. The
     * function returns the vectors containing the values of the functions
     * \f$f'(\chi)\f$ and \f$l(\chi)\f$ defined in Eqs. (30) and (31) of
     * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     * 
     * @param chi_vec Values of \f$\chi = k T\f$ at which the shape function
     *     is evaluated
     * @return Pair of vectors containing the values of \f$f'(\chi)\f$ and
     *     \f$l(\chi)\f$
     */
    std::pair<std::vector<double>,std::vector<double>> calculate_A_squared(
        const std::vector<double>& chi_vec
    );

    /**
     * @brief Calculate the shape function \f$|\mathcal{A}_+(\chi)|^2\f$,
     *     Eq. (35) of
     *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     * 
     * The shape function is evaluated for the values previously set using
     * @ref set_chi (or used in previous calls of @ref calculate_A_squared)
     * and can be retrieved calling @ref get_chi. The values of
     * \f$|\mathcal{A}_+(\chi)|^2\f$ are stored interval and can be retrieved
     * calling @ref get_shape_function. The function returns the vectors
     * containing the values of the functions \f$f'(\chi)\f$ and \f$l(\chi)\f$ 
     * defined in Eqs. (30) and (31) of
     * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     * 
     * @return Pair of vectors containing the values of \f$f'(\chi)\f$ and
     *     \f$l(\chi)\f$
     */
    std::pair<std::vector<double>,std::vector<double>> calculate_A_squared();

    /**
     * @brief Calculate the kinetic spectrum (momentum dependence of the
     *     stationary velocity field UETC).
     * 
     * The kinetic spectrum is evaluated on @p N logarithmically spaced values
     * between @p qR_min and @p qR_max. The values of \f$q R_*\f$ are stored
     * internally and can be retrieved calling @ref get_qR.
     *
     * @param qR_min Minimal value of \f$q R_*\f$
     * @param qR_max Maximal value of \f$q R_*\f$
     * @param N Number of points
     */
    void calculate_kinetic_spectrum(double qR_min, double qR_max, size_t N);

    /**
     * @brief Calculate the kinetic spectrum (momentum dependence of the
     *     stationary velocity field UETC).
     * 
     * The kinetic spectrum is evaluated for the values provided in @p qR_vec.
     * The values of \f$q R_*\f$ are stored internally and can be retrieved
     * calling @ref get_qR.
     * 
     * @param qR_vec Values of \f$q R_*\f$ at which the kinetic spectrum is
     *     evaluated
     */
    void calculate_kinetic_spectrum(const std::vector<double>& qR_vec);

    /**
     * @brief Calculate the kinetic spectrum (momentum dependence of the
     *     stationary velocity field UETC).
     * 
     * The kinetic spectrum is evaluated for the values previously set using
     * @ref set_qR (or used in previuos calls of
     * @ref calculate_kinetic_spectrum) and can be retrieved calling
     * @ref get_qR.
     */
    void calculate_kinetic_spectrum();

    /**
     * @brief Calculate the gravitational wave spectrum.
     * 
     * The gravitational wave spectrum is evaluated on @p N logarithmically
     * spaced values between @p kR_min and @p kR_max. The values of \f$k R_*\f$
     * are stored internally and can be retrieved calling @ref get_kR.
     * 
     * @param kR_min Minimal value of \f$k R_*\f$
     * @param kR_max Maximal value of \f$k R_*\f$
     * @param N Number of points
     * @param N_qR If provided, @ref set_qR_from_kR is called to set \f$q R_*\f$
     *     with @p N_qR points
     */
    void calculate_gw_spectrum(
        double kR_min, 
        double kR_max, 
        size_t N, 
        size_t N_qR=static_cast<size_t>(-1)
    );

    /**
     * @brief Calculate the gravitational wave spectrum.
     * 
     * The gravitational wave spectrum is evaluated for the values provided in
     * @p kR_vec. The values of \f$k R_*\f$ are stored internally and can be
     * retrieved calling @ref get_kR.
     * 
     * @param kR_vec Values of \f$k R_*\f$ at which the gravitational wave
     *     spectrum is evaluated
     * @param N_qR If provided, @ref set_qR_from_kR is called to set \f$q R_*\f$
     *     with @p N_qR points
     */
    void calculate_gw_spectrum(
        const std::vector<double>& kR_vec, size_t N_qR=static_cast<size_t>(-1)
    );

    /**
     * @brief Calculate the gravitational wave spectrum.
     * 
     * The gravitational wave spectrum is evaluated for the values previously
     * set using @ref set_kR (or used in previuos calls of 
     * @ref calculate_gw_spectrum) and can be retrieved calling @ref get_kR.
     */
    void calculate_gw_spectrum();

    /**
     * @brief Calculate the root-mean-squared fluid velocity squared
     *     \f$\bar{U}_f^2 = \langle u^2(x)\rangle\f$, cf. footnote 6 of
     *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     */
    void calculate_mean_squared_fluid_velocity();
    
    /** 
     * @brief Calculate the normalization \f$\mathcal{C}\f$ of the anisotropic
     *     stress UETC.
     */
    void calculate_anisotropic_stress_UETC_normalization();

    /**
     * @brief Calculate the anisotropic stress autorcorrelator (momentum
     *     dependence of the anisotropic stress UETC).
     * 
     * The anisotropic stress is evaluated at the same momenta as the
     * gravitational wave spectrum, which can be set and retrieved calling
     * @ref set_kR and @ref get_kR, respectively.
     */
    void calculate_anisotropic_stress_autocorrelator();    

    /**
     * @brief Calculate the spectral-function-weighted average of the GW
     *     amplitude dependence on the source duration.
     */
    void calculate_weighted_average_duration_dependence();

    /** @brief Get the @ref FluidProfile. */
    const FluidProfile& get_fluid_profile() const;

    /** 
     * @brief Get the \f$\chi = k T\f$ vector at which the shape function
     *     \f$|\mathcal{A}_+|^2\f$ is evaluated.
     */
    const std::vector<double>& get_chi() const;

    /**
     * @brief Get the \f$q R_*\f$ vector at which the kinetic spectrum is
     *     evaluated. 
     */
    const std::vector<double>& get_qR() const;

    /**
     * @brief Get the \f$k R_*\f$ vector at which the gravitational wave
     *     spectrum is evaluated. 
     */
    const std::vector<double>& get_kR() const;

    /**
     * @brief Get the shape function \f$|\mathcal{A}_+|^2\f$, Eq. (35) of
     *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)". 
     * 
     * The corresponding values of \f$\chi\f$ can be obtained with @ref get_chi.
     */
    const std::vector<double>& get_shape_function();

    /**
     * @brief Get the dimensionelss kinetic spectrum
     *     \f$R_*^{-1} E_\mathrm{kin}(q R_*)\f$ (i.e. the time-independent
     *     component of the velocity field UETC spectrum), cf. Eq. (39) of
     *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     * 
     * The corresponding values of \f$q R_*\f$ can be obtained with @ref get_qR.
     */
    const std::vector<double>& get_kinetic_spectrum();

    /**
     * @brief Get the normalized kinetic spectrum \f$\zeta_\mathrm{kin} =
     *     \frac{E_\mathrm{kin}}{E_\mathrm{kin}^\mathrm{peak}}\f$,
     *     defined in Eq. (41) of
     *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     * 
     * The corresponding values of \f$q R_*\f$ can be obtained with @ref get_qR.
     */
    std::vector<double> get_normalized_kinetic_spectrum();

    /** 
     * @brief Get the dimensionless spectral density of the plane wave
     *     components of the velocity field, \f$R_*^{-3} P_v(q R_*) =
     *     \frac{2 \pi^2}{R_*^2 q^2}\,R_*^{-1} E_\mathrm{kin}(q R_*) \f$, cf.
     *     Eq. (4.17) of @ref ref_Hindmarsh_2019 "Hindmarsh and Hijazi (2019)".
     * 
     * The corresponding values of \f$q R_*\f$ can be obtained with @ref get_qR.
     */
    std::vector<double> get_plane_wave_velocity_spectral_density();

    /** 
     * @brief Get the power spectrum, of the plane wave components of the
     *     velocity field, \f$\mathcal{P}_v(q) = \frac{q^3}{2 \pi^2} P_v(q)
     *     = q E_\mathrm{kin}(q)\f$.
     * 
     * The corresponding values of \f$q R_*\f$ can be obtained with @ref get_qR.
     */
    std::vector<double> get_plane_wave_velocity_power_spectrum();

    /** 
     * @brief Get the dimensionless spectral density of the velocity field,
     *     \f$R_*^{-3} P_\tilde{v}(q R_*) = 2 R_*^{-3} P_v(q R_*) =
     *     \frac{4 \pi^2}{(q R_*)^2}\,R^{-1} E_\mathrm{kin}(q R_*)\f$.
     * 
     * The corresponding values of \f$q R_*\f$ can be obtained with @ref get_qR.
     */
    std::vector<double> get_velocity_spectral_density();

    /** 
     * @brief Get the velocity power spectrum \f$\mathcal{P}_{\tilde{v}}(q R_*)
     *     = \frac{q^3}{2 \pi^2} P_\tilde{v}(q) = 2 q E_\mathrm{kin}(q)\f$, cf.
     *     Eq. (4.18) of @ref ref_Hindmarsh_2019 "Hindmarsh and Hijazi (2019)".
     * 
     * The corresponding values of \f$q R_*\f$ can be obtained with @ref get_qR.
     */
    std::vector<double> get_velocity_power_spectrum();

    /** 
     * @brief Get the gravitational wave spectrum \f$\frac{d\Omega_\mathrm{GW}^*
     *     }{d \log k}\f$ at production, cf. Eq. (15) and (93) of
     *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)"
     *     (with \f$\mathcal{T}_\mathrm{GW}=1\f$).
     * 
     * The corresponding values of \f$k R_*\f$ can be obtained with @ref get_kR.
     */
    const std::vector<double>& get_gravitational_wave_spectrum();

    /** 
     * @brief Get the gravitational wave power spectrum
     *     \f$\mathcal{P}_\mathrm{gw}(k)\f$, cf. Eq. (3.6) of
     *     @ref ref_Hindmarsh_2019 "Hindmarsh and Hijazi (2019)".
     * 
     * This is the same as @ref get_gravitational_wave_spectrum.
     * 
     * The corresponding values of \f$k R_*\f$ can be obtained with @ref get_kR.
     */
    const std::vector<double>& get_gravitational_wave_power_spectrum();

    /** 
     * @brief Get the dimensionless gravitational wave spectral density
     *     \f$R_*^{-3} P_\mathrm{gw}(k R_*)\f$.
     * 
     * The corresponding values of \f$k R_*\f$ can be obtained with @ref get_kR.
     */
    std::vector<double> get_gravitational_wave_spectral_density();

    /** 
     * @brief Get the anisotropic stress autocorrelator power spectrum
     *     \f$\mathcal{P}_\Pi(k) = k E_\Pi(k)\f$, cf. Eq. (44) of
     *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     * 
     * The corresponding values of \f$k R_*\f$ can be obtained with @ref get_kR.
     */
    const std::vector<double>& get_anisotropic_stress_power_spectrum();

    /** 
     * @brief Get the normalized anisotropic stress autocorrelator power
     *     spectrum \f$(k R_*)^3 \zeta_\Pi(k)\f$, cf. Eq. (47) of
     *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     * 
     * The corresponding values of \f$k R_*\f$ can be obtained with @ref get_kR.
     */
    std::vector<double> get_normalized_anisotropic_stress_power_spectrum();

    /** 
     * @brief Get the dimensionless anisotropic stress autocorrelator
     *     spectral density
     *     \f$R_*^{-3} P_\Pi(k) = \frac{2\pi^2}{(R_* k)^2} R_*^{-1} E_\Pi(k)\f$.
     * 
     * The corresponding values of \f$k R_*\f$ can be obtained with @ref get_kR.
     */
    std::vector<double> get_anisotropic_stress_spectral_density();

    /** 
     * @brief Get the normalized gravitational wave spectrum
     *     \f$(k R_*)^3 \zeta_\mathrm{GW}(k R_*)\f$, cf. Eq. (93) of
     *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     * 
     * The corresponding values of \f$k R_*\f$ can be obtained with @ref get_kR.
     */
    std::vector<double> get_normalized_gravitational_wave_spectrum();

    /**
     * @brief Get the spectral modification function \f$\tilde{\Delta} =
     *     \frac{\zeta_\mathrm{GW}}{\zeta_\Pi}\f$, cf. Eq. (95) of
     *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     * 
     * The corresponding values of \f$k R_*\f$ can be obtained with @ref get_kR.
     */
    std::vector<double> get_spectral_modification_function();

    /**
     * @brief Get the root-mean-squared fluid velocity squared
     *     \f$\bar{U}_f^2 = \langle u^2(x)\rangle\f$, cf. footnote 6 of
     *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     */
    double get_mean_squared_fluid_velocity();

    /** 
     * @brief Get the total kinetic energy density parameter
     *     \f$\Omega_K = \Gamma \bar{U}_f^2\f$ of the kinetic energy density,
     *     where \f$\Gamma\f$ is the adiabatic index of the fluid. 
     */
    double get_kinetic_energy_density_parameter();

    /**
     * @brief Get the peak amplitude \f$R_*^{-1} E_\mathrm{kin}^*\f$ of the
     *     kinetic spectrum, cf. Eq. (41) of
     *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     */
    double get_kinetic_spectrum_amplitude();

    /**
     * @brief Get the peak position \f$(k R_*)_\mathrm{kin}^\mathrm{peak}\f$
     *     of the kinetic spectrum. 
     */
    double get_kinetic_spectrum_peak_position();

    /** 
     * @brief Get the normalized kinetic energy density parameter
     *     \f$\mathcal{K} = \frac{R_* \bar{U}_f^2}{2\,E_\mathrm{kin}^*}\f$,
     *     cf. Eqs. (42) and (43) of
     *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     */
    double get_normalized_kinetic_density_parameter();

    /**
     * @brief Get the normalization \f$\mathcal{C}\f$ of the anisotropic
     *     stress spectral shape, cf. Eq. (46) of
     *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     */
    double get_anisotropic_stress_UETC_normalization();

    /**
     * @brief Get the peak amplitude \f$\left(K^3 \zeta_\Pi\right)_\mathrm{peak}
     *     \f$ of the normalized anisotropic stress spectrum.
     */
    double get_normalized_anisotropic_stress_UETC_amplitude();

    /**
     * @brief Get the peak position \f$(k R_*)_\mathrm{\Pi}^\mathrm{peak}\f$ of
     *     the normalized anisotropic stress spectrum. 
     */
    double get_normalized_anisotropic_stress_UETC_peak_position();

    /**
     * @brief Get the spectral-function-weighted average
     *     \f$\tilde{\Delta}_0(\delta\tau_\mathrm{fin}, R_*)\f$ of GW amplitude 
     *     dependence on the source duration, Eq. (63) of
     *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
     */
    double get_weighted_average_duration_dependence();

    /**
     * @brief Set the nucleation history type.
     * @param nhist Nucleation history (exponential or simultaneous)
     */
    void set_nucleation_history(sound_shell_model::NucleationHistory nhist);

    /** @brief Get the nucleation history type.  */
    sound_shell_model::NucleationHistory get_nucleation_history() const;

    /** 
     * @brief Change the wall velocity and/or stength parameter and recalculate
     *     the fluid profile. 
     *
     * @param v_w wall veclocity \f$\xi_w\f$
     * @param alpha_N strength parameter \f$\alpha_N\f$ (far in front of the
     *     wall)
     * @param step_size Step size for the fluid velocity (or coordinate) when
     *     solving the differential equation (default: 1e-4)
     * @param solve_mode \c profile_calculator::ProfileSolveMode switch
     *     determining which differential equation for the velocity profile is
     *     solved (default: 
     *     \c profile_calculator::ProfileSolveMode::D_V_D_XI_FOR_SHOCKS)
     */
    void recalculate_profile(
        double v_w, 
        double alpha_N, 
        double step_size=1e-3,
        profile_calculator::ProfileSolveMode solve_mode
            = profile_calculator::ProfileSolveMode::D_V_D_XI_FOR_SHOCKS
    );

    /** @brief Get the value of the wall velocity \f$\xi_w\f$. */
    double get_wall_velocity() const;

    /** @brief Get the value of the strength parameter \f$\alpha_N\f$. */
    double get_transition_strength() const;

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

    /** 
     * @brief Set the expansion history type. 
     * @param ehist Expansion history (radiation domination or static, see
     *     @ref sound_shell_model::ExpansionHistory)
     */
    void set_expansion_history(sound_shell_model::ExpansionHistory ehist);

    /** @brief Get the expansion history type. */
    sound_shell_model::ExpansionHistory get_expansion_history() const;

    /**
     * @brief Set the mean bubble separation \f$H_* R_*\f$.
     * @param R_H Mean bubble separation at collision times Hubble rate,
     *     \f$H_* R_*\f$
     */
    void set_mean_bubble_separation(double R_H);

    /** 
     * @brief Get the mean bubble separation at collision times Hubble rate,
     *     \f$H_* R_*\f$
     */
    double get_mean_bubble_separation() const;

    /** 
     * @brief Set the source duration to a fixed value.
     * @param duration Source duration \f$H_* \delta\tau_\mathrm{fin}\f$ in
     *     units of the Hubble time
     */
    void set_source_duration(double duration);

    /** 
     * @brief Calulcate the source duration.
     * 
     * This function calculates the source duration from the mean bubble
     * separation \f$H_* R_*\f$ and the kinetic energy density \f$\Omega_K\f$,
     * \f$H_* \delta\tau_\mathrm{fin} = \frac{H_* R_*}{\sqrt{\Omega_K}}\f$.
     */
    void calculate_source_duration();

    /** 
     * @brief Get the duration of the sound wave source in units of the Hubble
     *     time.
     */
    double get_source_duration();

    /** 
     * @brief Set settings for the momentum integration of anisotropic stress
     *     correlators.
     * @param n_p Number of points in the integration over \f$p R_*\f$
     * @param n_z Number of points in the integration over \f$z\f$ 
     * @param integ Approximation used in the momentum integration, cf. @ref
     * sound_shell_model::AnisotropicStressCorrelatorSettings::MomentumIntegration
     */
    void set_correlator_integration_settings(
        size_t n_p, 
        size_t n_z, 
        sound_shell_model::MomentumIntegration integ
            = sound_shell_model::MomentumIntegration::FULL
    );

private:
    /// verbosity flag
    inline static bool verbose = false; 
    /// fluid profile from which the spectra are calculated
    FluidProfile fluid_profile; 
    /// nucleation history
    sound_shell_model::NucleationHistory nuc_hist; 
    /// expansion history
    sound_shell_model::ExpansionHistory exp_hist;  
    /// whether the source duration is calculated or set manually
    bool calculate_duration = true;
    /// source duration times Hubble rate \f$H_* \delta \tau_\mathrm{fin}\f$
    double H_delta_tau_fin;
    /// bubble separation at collision times Hubble rate, \f$H_* R_*\f$
    double H_R;
    /// root-mean-squared fluid velocity squared \f$\bar{U}_f^2\f$
    double Uf2 = -1.0;
    /**
     * @brief normalization \f$\mathcal{C}\f$ of the anisotropic stress
     * autocorrelator power spectrum
     */
    double script_C = -1.0;
    /**
     * @brief Kinetic-spectrum-weighted average of the duration dependence of
     *     the GW spectrum
     */
    double Delta_tilde_0 = -1.0;
    /// dimensionless momentum \f$k R_*\f$ at which the GW spectrum is evaluated
    std::vector<double> kR;
    /**
     * @brief dimensionless momentum \f$q R_*\f$ at which the kinetic spectrum
     *     is evaluated
     */
    std::vector<double> qR;
    /**
     * @brief dimensionless momentum \f$\chi = k T\f$ at which the sine
     *     transformed profiles are evaluated
     */
    std::vector<double> chi;
    /// shape function \f$|\mathcal{A}_+|^2\f$
    std::vector<double> A2;
    /**
     * @brief dimensionless kinetic spectrum \f$R_*^{-1} E_\mathrm{kin}(q)\f$
     *     (momentum dependence of the stationary velocity field UETC)
     */
    std::vector<double> Ekin_over_R;
    /// Anisotropic stress autocorrelator \f$k E_\Pi\f$
    std::vector<double> k_EPi;
    /// GW spectrum \f$\Omega_\mathrm{GW}^*\f$ at production
    std::vector<double> OmegaGW;
    /**
     * @brief Number of points in the integration over \f$p R_*\f$ when
     *     calculating anisotropic stress correlators
     */
    size_t N_p = 5000;
    /**
     * @brief Number of points in the integration over \f$z\f$ when calculating
     *     anisotropic stress correlators
     */
    size_t N_z = 100;
    /**
     * @brief Approximation used in the momentum integration of anisotropic
     *     stress correlators, cf. @ref
     * sound_shell_model::AnisotropicStressCorrelatorSettings::MomentumIntegration
     */
    sound_shell_model::MomentumIntegration integration
        = sound_shell_model::MomentumIntegration::FULL;

}; // class SoundShellSpectrum

} // namespace inverse_pt