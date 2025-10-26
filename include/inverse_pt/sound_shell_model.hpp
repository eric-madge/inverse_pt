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

#include "inverse_pt/utils.hpp"
#include "inverse_pt/constants.hpp"

#include <cmath>
#include <tuple>
#include <vector>

/**
 * @namespace sound_shell_model
 * @brief Functions to calculate the gravitational wave spectrum within the 
 *     sound shell model.
 * 
 * This is based on
 * 1. @anchor ref_Hindmarsh_2018
 *    M. Hindmarsh, 
 *    "Sound shell model for acoustic gravitational wave production at a
 *    first-order phase transition in the early Universe," 
 *    [Phys. Rev. Lett. **120** (2018) no. 7, 071301][1], 
 *    [arXiv: 1608.04735 [astro-ph.CO]](https://arxiv.org/abs/1608.04735).
 * 2. @anchor ref_Hindmarsh_2019
 *    M. Hindmarsh and M. Hijazi,
 *    "Gravitational waves from first order cosmological phase transitions in
 *    the Sound Shell Model,"
 *    [JCAP **12** (2019), 062](https://doi.org/10.1088/1475-7516/2019/12/062),
 *    [arXiv:1909.10040 [astro-ph.CO]](https://arxiv.org/abs/1909.10040).
 * 3. @anchor ref_Roper_Pol_2024
 *    A. Roper Pol, S. Procacci and C. Caprini,
 *    "Characterization of the gravitational wave spectrum from sound waves
 *    within the sound shell model,"
 *    [Phys. Rev. D **109** (2024) 6, 063531][2],
 *    [arXiv:2308.12943 [gr-qc]](https://arxiv.org/abs/2308.12943).
 * 
 * [1]: https://doi.org/10.1103/PhysRevLett.120.071301
 * [2]: https://doi.org/10.1103/PhysRevD.109.063531
 */
namespace inverse_pt::sound_shell_model {

/// Nucleation history options
enum class NucleationHistory{
    EXPONENTIAL,  ///< exponential nucleation
    SIMULTANEOUS  ///< simultaneous nucleation
};

/// Expansion History options
enum class ExpansionHistory{
    RADIATION_DOMINATION, ///< assume a radiation dominated universe
    STATIC_UNIVERSE,      ///< assume a static (non-expanding) universe
    /// no expansion history / time integration (i.e. use \f$\Delta = 1\f$ to 
    /// calculate the autocorrelation \f$\zeta_\Pi\f$)
    NONE
};

/**
 * @brief Settings for the momentum and time integration of the anisotropic
 * stress and gravitational wave spectrum in
 * @ref anisotropic_stress_correlator_integration.
 */
struct AnisotropicStressCorrelatorSettings {
    /// Approximation used in the momentum integration
    enum class MomentumIntegration {
        /** 
         * Full momentum integration, using Eqs. (56) and (59) of 
         * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)"
         * for \f$\Delta\f$ and \f$\Delta^\mathrm{flat}\f$, respectively.
         */
        FULL,
        /**
         * Use the low-frequency limit
         * \f$\Delta_0\f$ or \f$\Delta_0^\mathrm{flat}\f$, i.e. Eqs. (61,62) of
         * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)" for
         * \f$\Delta\f$ and \f$\Delta^\mathrm{flat}\f$.
         */
        LOW_FREQUENCY,
        /**
         * Use the inifinite-duration (high-frequency) limit of \f$\Delta\f$,
         * \f$k\delta\tau_\mathrm{fin}\gg 1\f$, i.e. Eqs. (86) and (83) of
         * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)",
         * respectively. In this limit, \f$\Delta\f$ reduces to Dirac
         * \f$\delta\f$-functions, simplifying the momentum integration.
         */
        INFINITE_DURATION
    };

    /// Number of points in the integration over \f$p R_*\f$ (default: 5000)
    size_t N_p = 5000;
    /// Number of points in the integration over \f$z\f$ (default: 100) 
    size_t N_z = 100;  
    /** 
     * @brief Source duration in units of the Hubble time,
     *    \f$\mathcal{H}_* \delta\tau_\mathrm{fin}\f$
     */
    double H_delta_tau = 0.0;
    /**
     * @brief Mean bubble separation in units of the Hubble size,
     *     \f$\mathcal{H}_* R_*\f$
     */
    double H_R = 1.0; 
    /// Assumptions used in the time integration 
    ExpansionHistory universe = ExpansionHistory::NONE;
    /// Approximation used in the momentum integration of the time integrals
    MomentumIntegration freq_approx = MomentumIntegration::FULL;
    /** 
     * @brief Numeric prefactor of the integral in 
     *     @ref anisotropic_stress_correlator_integration (\f$2 \bar{w}^2\f$
     *     for \f$k E_\Pi\f$ and \f$3\Gamma^2\f$ for \f$\Omega_\mathrm{GW}\f$)
     */
    double c_w = 2.0;

    /**
     * @brief Default constructor (for calculating the anisotropic stress
     *     autocorrelation power spectrum \f$k E_\Pi(\tau,\tau,k)\f$.
     */
    AnisotropicStressCorrelatorSettings() = default;
    /** 
     * @brief Constuctor for settings for calculating the anisotropic stress
     *     autocorrelation power spectrum \f$k E_\Pi(\tau,\tau,k)\f$.
     * @note This sets @ref universe to \c ExpansionHistory::None and 
     *     @ref freq_approx to \c MomentumIntegration::FULL.
     * @param N_p_ Number of points in the integration over \f$p R_*\f$
     * @param N_z_ Number of points in the integration over \f$z\f$ 
     * @param c_w_ Numeric prefactor of the integral in
     *     @ref anisotropic_stress_correlator_integration
     *     (optional, default: 2.0)
     */
    AnisotropicStressCorrelatorSettings(
        size_t N_p_, size_t N_z_, double c_w_=2.0
    ) : 
        N_p(N_p_), 
        N_z(N_z_), 
        universe(ExpansionHistory::NONE), 
        freq_approx(MomentumIntegration::FULL), 
        c_w(c_w_)
    {};
    /** 
     * @brief Constructor for time integration settings for the graviational
     *     wave spectrum \f$\Omega_\mathrm{GW}\f$.
     * @note This automatically sets @ref c_w to `16.0/3.0` (i.e.
     *     \f$3 \Gamma^2\f$ for \f$\Gamma = \frac{\bar{w}}{\bar{e}} =
     *     \frac{4}{3}\f$).
     * @param H_delta_tau_ Source duration in units of the Hubble time,
     *     \f$\mathcal{H}_* \delta\tau_\mathrm{fin}\f$
     * @param H_R_ Mean bubble separation in units of the Hubble size,
     *     \f$\mathcal{H}_* R_*\f$
     * @param universe_ Assumptions used in the time integration, cf.
     *     @ref ExpansionHistory (optional, default:
     *     \c ExpansionHistory::RADIATION_DOMINATION)
     * @param freq_approx_ Approximation used in the momentum integration of the
     *     time integrals, cf. 
     *     @ref AnisotropicStressCorrelatorSettings::MomentumIntegration
     *     (optional, default: \c MomentumIntegration::FULL,
     *     i.e. no approximation)
     * @param c_w_ Numeric prefactor of the integral in
     *     @ref anisotropic_stress_correlator_integration
     *     (optional, default: 16.0/3.0)
     */
    AnisotropicStressCorrelatorSettings(
        double H_delta_tau_, 
        double H_R_, 
        ExpansionHistory universe_=ExpansionHistory::RADIATION_DOMINATION, 
        MomentumIntegration freq_approx_= MomentumIntegration::FULL,
        double c_w_=16.0/3.0
    ) : 
        H_delta_tau(H_delta_tau_), 
        H_R(H_R_), 
        universe(universe_), 
        freq_approx(freq_approx_),
        c_w(c_w_)
    {};
    /** 
     * @brief Full constructor.
     * @param N_p_ Number of points in the integration over \f$p R_*\f$
     * @param N_z_ Number of points in the integration over \f$z\f$ 
     * @param H_delta_tau_ Source duration in units of the Hubble time,
     *     \f$\mathcal{H}_* \delta\tau_\mathrm{fin}\f$
     * @param H_R_ Mean bubble separation in units of the Hubble size,
     *     \f$\mathcal{H}_* R_*\f$
     * @param universe_ Assumptions used in the time integration, cf.
     *     @ref ExpansionHistory (optional, default:
     *     ExpansionHistory::RADIATION_DOMINATION)
     * @param freq_approx_ Approximation used in the momentum integration of the
     *     time integrals, cf. 
     *     @ref AnisotropicStressCorrelatorSettings::MomentumIntegration
     *     (optional, default: \c MomentumIntegration::Full,
     *     i.e. no approximation)
     * @param c_w_ Numeric prefactor of the integral in 
     *     @ref anisotropic_stress_correlator_integration
     *     (optional, default: 16.0/3.0)
     */
    AnisotropicStressCorrelatorSettings(
        size_t N_p_, 
        size_t N_z_, 
        double H_delta_tau_, 
        double H_R_, 
        ExpansionHistory universe_=ExpansionHistory::RADIATION_DOMINATION, 
        MomentumIntegration freq_approx_= MomentumIntegration::FULL,
        double c_w_=16.0/3.0
    ) : 
        N_p(N_p_), 
        N_z(N_z_), 
        H_delta_tau(H_delta_tau_), 
        H_R(H_R_), 
        universe(universe_), 
        freq_approx(freq_approx_), 
        c_w(c_w_) 
    {};
};

using MomentumIntegration = 
    AnisotropicStressCorrelatorSettings::MomentumIntegration;

/**
 * @brief Bubble lifetime distribution function for exponential nucleation.
 * 
 * This is Eq. (4.27) of @ref ref_Hindmarsh_2019 "Hindmarsh and Hijazi (2019)".
 * 
 * @param T Dimensionless bubble lifetime \f$\tilde{T} = \beta T_i\f$
 * @return Bubble lifetime distribution evaluated at T for exponential
 *     nucleation
 */
#pragma omp declare simd
inline double nu_exp(double T) { return std::exp(-T); }

/**
 * @brief Bubble lifetime distribution function for simultaneous nucleation.
 * 
 * This is Eq. (4.32) of @ref ref_Hindmarsh_2019 "Hindmarsh and Hijazi (2019)".
 * 
 * @param T Dimensionless bubble lifetime \f$\tilde{T} = \beta T_i\f$
 * @return Bubble lifetime distribution evaluated at T for simultaneous
 *     nucleation
 */
#pragma omp declare simd
inline double nu_sim(double T) { 
    return 0.5 * T * T * std::exp(- T * T * T / 6.0); 
}

/** 
 * @brief Calculate the (derivative of) sine transforms of the velocity and
 *     energy spectrum.
 * 
 * This function calculates the derivative of the sine transform \f$f(\chi)\f$
 * of the velocity spectrum \f$v\f$, the sine transform \f$l(\chi)\f$ of the
 * energy fluctuation variable \f$\lambda = \frac{e-\bar{e}}{\bar{w}}\f$, and
 * the absolute value squared of their sum, \f$|\mathcal{A}_+(\chi)|^2 
 * = \frac{1}{4} \left[(f'(\chi))^2 + (c_s l(\chi))^2\right]\f$, given in Eqs.
 * (29-31) of @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)"
 * [or Eqs. (4.5), (4.11) and (4.8) of
 * @ref ref_Hindmarsh_2019 "Hindmarsh and Hijazi (2019)" with
 * \f$A=\mathcal{A}_+\f$],
 * \f[
 *     f'(\chi) = \frac{d}{d \chi} \left[ \frac{4 \pi}{\chi} 
 *     \int\limits_0^\infty d\xi\,v(\xi) \sin(\chi \xi) \right]
 *     \qquad \text{and} \qquad
 *     l(\chi) = \frac{4 \pi}{\chi} \int\limits_0^\infty d\xi\,
 *     \frac{e-\bar{e}}{\bar{w}} \xi \sin(\chi \xi) \,,
 * \f]
 * where \f$\bar{e}\f$ and \f$\bar{w}\f$ are the mean energy and enthalpy
 * density, and \f$\chi = k T_n\f$ with \f$T_n\f$ being the lifetime of the 
 * \f$n\f$-th bubble.
 * 
 * @param chi_vec Array of values of the dimensionless wavenumber
 *     \f$\chi = k T_n\f$
 * @param xi Values of the self-similar coordinate \f$\xi\f$
 * @param v Fluid velocity profile \f$v(\xi)\f$
 * @param w Enthalpy density profile \f$w(\xi)\f$
 * @param e Energy density profile \f$e(\xi)\f$
 * @param idx_wall Index of the bubble wall position
 *     (first element outside the wall)
 * @param idx_sh Index of the shock front position (first element outside the
 *     shock, use \c constants::NO_INDEX if there is no shock front)
 * @return The sine transformed functions \f$|A(\chi)|^2\f$, \f$f'(\chi)\f$ and
 *     \f$l(\chi)\f$
*/
auto calculate_sine_transforms(
    const std::vector<double>& chi_vec, 
    const std::vector<double>& xi, 
    const std::vector<double>& v, 
    const std::vector<double>& w, 
    const std::vector<double>& e, 
    size_t idx_wall, 
    size_t idx_sh=constants::NO_INDEX
) -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>;

/**
 * @brief Scaled kinetic spectrum \f$\bar{E}_\mathrm{kin}(k R_*) = R_*^{-1} 
 *     E_\mathrm{kin}(k)\f$.
 * 
 * This is the momentum dependent coefficient \f$E_\mathrm{kin}(k) = 
 * E_\mathrm{kin}^{(1)}(k)\f$ of the velocity field unequal-time correlator,
 * \f[
 *     E_\mathrm{kin}(k, \tau_1, \tau_2) =
 *      E_\mathrm{kin}(k) \cos\left(k c_s (tau_1-\tau_2)\right),
 * \f]
 * cf. Eqs. (39), (33) and (35) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)". We use the
 * dimensionless momentum \f$z = k R_*\f$ as input variable and divide by $R_*$
 * to make the kinetic spectrum dimensionless.
 * 
 * Note that this is related to the velocity spectral density of the plane wave
 * components of the velocity field, Eq. (4.17) of
 * @ref ref_Hindmarsh_2019 "Hindmarsh and Hijazi (2019)",
 * \f$P_v(k) = \frac{2 \pi^2}{k^2} E_\mathrm{kin}(k)\f$.
 * 
 * @param kR Momentum \f$k\f$ times bubble separation \f$R_*\f$ at which the
 *     kinetic spectrum is evaluated
 * @param chi_vec Values of \f$\chi = T q = \tilde{T} q / \beta\f$ at which
 *     \f$|\mathcal{A}_+|^2\f$ is provided
 * @param A2_vec Precomputed values of the shape function
 *     \f$|\mathcal{A}_+|^2\f$ used in the integration
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param nuc_hist Nucleation history (exponential of simultaneous) used for the
 *     collision time distribution (default: exponential nucleation)
 * @return Scaled kinetic spectrum \f$R_*^{-1} E_\mathcal{kin}(k R_*)\f$
 */
std::vector<double> scaled_kinetic_spectrum(
    const std::vector<double>& kR, 
    const std::vector<double>& chi_vec, 
    const std::vector<double>& A2_vec, 
    double xi_wall, 
    NucleationHistory nuc_hist=NucleationHistory::EXPONENTIAL
);

/**
 * @brief Spectral shape \f$\zeta_\mathrm{kin}(k R_*)\f$ of the kinetic
 *     spectrum.
 * 
 * This is the spectral shape \f$\zeta_\mathrm{kin}(k R_*)\f$ of the kinetic
 * spectrum, \f$E_\mathrm{kin}(k) = E_\mathrm{kin}^*
 * \zeta_\mathrm{kin}(k R_*)\f$ (cf. Eq. (41) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)"), where
 * \f$E_\mathrm{kin}^*\f$ is the maximum amplitude of the kinetic spectrum and
 * \f$R_*\f$ is the mean bubble separation.
 * 
 * @param E_kin_over_R Scaled kinetic spectrum
 *     \f$\bar{E}_\mathrm{kin}(k R_*) = R_*^{-1} E_\mathrm{kin}(k)\f$
 * @return kinetic spectral shape \f$\zeta_\mathrm{kin}(k R_*) =
 *     \frac{E_\mathrm{kin}(k)}{E_\mathrm{kin}^*}\f$
 */
std::vector<double> kinetic_spectral_shape(
    const std::vector<double>& E_kin_over_R
);

/**
 * @brief Mean square fluid velocity squared \f$\bar{U}_f^2\f$.
 * 
 * This is the root-mean-square fluid velocity squared, cf. Eqs. (4.33) and
 * (4.34) of @ref ref_Hindmarsh_2019 "Hindmarsh and Hijazi (2019)", and footnote
 * 6 of @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)", 
 * calculated from the kinetic spectrum (or velocity-field UETC at equal times,
 * see also Eq. (42) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)")
 * \f[
 *     \bar{U}_f^2 = \langle \mathbf{u}^2(\mathbf{x}) \rangle
 *     = 2 \int\limits_{0}^{\infty} dk\, E_\mathrm{kin}(\tau, \tau, k)
 *     = 2\, \frac{E_\mathrm{kin}^*}{R_*} \mathcal{K} \,,
 * \f]
 * where \f$E_\mathrm{kin}^*\f$ is the maximum amplitude of the kinetic spectrum
 * and \f$R_*\f$ is the mean bubble separation.
 * 
 * @note The sound shell model assumes non-relativistic fluid velocities and a
 *     constant entropy density, \f$w(x) = \bar{w} \approx w_N\f$, where
 *     \f$w_N\f$ is the enthalphy density at the nucleation temperature. Thus,
 *     the root-mean-square fluid velocity coincides with the enthalpy-weighted
 *     average fluid velocity.
 * 
 * @param kR Momentum \f$k\f$ times bubble separation \f$R_*\f$ at which the
 *     kinetic spectrum is evaluated
 * @param E_kin_over_R Scaled kinetic spectrum
 *     \f$\bar{E}_\mathrm{kin}(k R_*) = R_*^{-1} E_\mathrm{kin}(k)\f$
 * @return Mean-square fluid velocity \f$\bar{U}_f^2\f$
 */
double mean_squared_fluid_velocity(
    const std::vector<double>& kR, const std::vector<double>& E_kin_over_R
);

/**
 * @brief Mean square fluid velocity squared \f$\bar{U}_f^2\f$.
 * 
 * This is the root-mean-square fluid velocity squared, cf. Eqs. (4.33) and
 * (4.34) of @ref ref_Hindmarsh_2019 "Hindmarsh and Hijazi (2019)". 
 * The integral over \f$\tilde{T}\f$ in Eq. (4.33) has been carried out
 * analytically, assuming that the nucleation history is either exponential or
 * simultaneous.
 * \f[
 *     \bar{U}_f^2 = \frac{3}{4 \pi \xi_w^3} \int\limits_0^\infty \!d\chi\,
 *     \frac{\chi^2}{2 \pi^2} \left|\mathcal{A}_+(\chi)\right|^2 \,.
 * \f]
 * 
 * @note The sound shell model assumes non-relativistic fluid velocities and a
 *     constant entropy density, \f$w(x) = \bar{w} \approx w_N\f$, where
 *     \f$w_N\f$ is the enthalphy density at the nucleation temperature. Thus,
 *     the root-mean-square fluid velocity coincides with the enthalpy-weighted
 *     average fluid velocity.
 * 
 * @param chi_vec Values of \f$\chi = T q = \tilde{T} q / \beta\f$ at which
 *     \f$|\mathcal{A}_+|^2\f$ is provided
 * @param A2_vec Precomputed values of the shape function
 *     \f$|\mathcal{A}_+|^2\f$ used in the integration
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @return Mean-square fluid velocity \f$\bar{U}_f^2\f$
 */
double mean_squared_fluid_velocity(
    const std::vector<double>& chi_vec, 
    const std::vector<double>& A2_vec,
    double xi_wall
);

/**
 * @brief Total kinetic energy density parameter \f$\Omega_K\f$.
 * 
 * This is the total kinetic energy density
 * \f$\Omega_K = \frac{\rho_\mathrm{kin}}{\rho_c}\f$
 * as a fraction of the critical energy density, where the kinetic energy is
 * given by \f$\rho_\mathrm{kin} = \langle w \gamma^2 v^2\rangle\f$. It is
 * related to the root-mean-square fluid velocity by 
 * \f$\Omega_K = \Gamma \bar{U}_f^2\f$, where \f$\Gamma = \frac{\bar{w}}{
 * \bar{e}}\f$ is the mean adiabatic index of the fluid in the stable phase.
 * 
 * @note This assumes \f$\Gamma=\frac{4}{3}\f$ by default. To obtain
 *     \f$\Omega_K = \frac{1}{2}\,\bar{U}_f^2\f$ as defined in Eq. (42) of 
 *     @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)", set
 *     @p Gamma to \c 0.5.
 * 
 * @param Uf2 mean-square fluid velocity \f$\bar{U}_f^2\f$
 * @param Gamma Mean adiabatic index \f$\Gamma = \frac{\bar{w}}{\bar{e}}\f$ of
 *     the fluid in the stable phase (default: 4/3)
 * @return kinetic energy density parameter \f$\Omega_K\f$
 */
double kinetic_energy_density_parameter(double Uf2, double Gamma=4.0/3.0);

/**
 * @brief Total kinetic energy density parameter \f$\Omega_K\f$.
 * 
 * This is the total kinetic energy density
 * \f$\Omega_K = \frac{\rho_\mathrm{kin}}{\rho_c}\f$
 * as a fraction of the critical energy density, calculated from the kinetic
 * spectrum (or velocity-field UETC at equal times, cf. Eq. (42) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)")
 * \f[
 *     \Omega_K = 2 \Gamma\, \int\limits_{0}^{\infty} dk\,
 *         E_\mathrm{kin}(\tau, \tau, k)
 *     = 2 \Gamma\, \frac{E_\mathrm{kin}^*}{R_*}\,\mathcal{K}\,,
 * \f]
 * where \f$Gamma=\frac{\bar{w}}{\bar{e}}\f$ is the mean adiabatic index of the 
 * fluid, \f$E_\mathrm{kin}^*\f$ is the maximum amplitude of the kinetic
 * spectrum and \f$R_*\f$ is the mean bubble separation. 
 * 
 * @note Set @p Gamma to `0.5` to obtain \f$\Omega_K\f$ as defined in Eq. (42)
 *     of @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
 * @note To obtain the integral of the kinetic spectral shape defined in Eq.
 *     (43) of @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)", 
 *     \f$\mathcal{K} = \int d K\,\zeta_\mathrm{kin}(K)\f$, provide the kinetic
 *     spectral shape \f$\zeta_\mathrm{kin}(k R_*)\f$ for @p E_kin_over_R and
 *     set @p Gamma to `0.5`.
 * 
 * @param kR Momentum \f$k\f$ times bubble separation \f$R_*\f$ at which the
 *     kinetic spectrum is evaluated
 * @param E_kin_over_R Scaled kinetic spectrum
 *     \f$\bar{E}_\mathrm{kin}(k R_*) = R_*^{-1} E_\mathrm{kin}(k)\f$
 * @param Gamma Mean adiabatic index \f$\Gamma = \frac{\bar{w}}{\bar{e}}\f$ of
 *     the fluid in the stable phase (default: 4/3)
 * @return kinetic energy density parameter \f$\Omega_K\f$
 */
double kinetic_energy_density_parameter(
    const std::vector<double>& kR, 
    const std::vector<double>& E_kin_over_R,
    double Gamma=4.0/3.0
);

/**
 * @brief Total kinetic energy density parameter \f$\Omega_K\f$.
 * 
 * This is the total kinetic energy density
 * \f$\Omega_K = \frac{\rho_\mathrm{kin}}{\rho_c}\f$
 * as a fraction of the critical energy density.
 * The integral over \f$\tilde{T}\f$ in Eq. (4.33) has been carried out
 * analytically, assuming that the nucleation history is either exponential or
 * simultaneous, as in the step from Eq. (4.33) to (4.34) in
 * @ref ref_Hindmarsh_2019 "Hindmarsh and Hijazi (2019)".
 * \f[
 *     \Omega_K = \frac{3}{4 \pi \xi_w^3}\, \Gamma\,
 *     \int\limits_0^\infty \!d\chi\,
 *     \frac{\chi^2}{2 \pi^2} \left|\mathcal{A}_+(\chi)\right|^2 \,.
 * \f]
 * 
 * @param xi_wall Wall velocity \f$\xi_w\f$
 * @param chi_vec Values of \f$\chi = T q = \tilde{T} q / \beta\f$ at which
 *     \f$|\mathcal{A}_+|^2\f$ is provided
 * @param A2_vec Precomputed values of the shape function
 *     \f$|\mathcal{A}_+|^2\f$ used in the integration
 * @param Gamma Mean adiabatic index \f$\Gamma = \frac{\bar{w}}{\bar{e}}\f$ of
 *     the fluid in the stable phase (default: 4/3)
 * @return kinetic energy density parameter \f$\Omega_K\f$
 */
double kinetic_energy_density_parameter(
    double xi_wall,
    const std::vector<double>& chi_vec,
    const std::vector<double>& A2_vec,
    double Gamma=4.0/3.0
);

/**
 * @brief Normalized amplitude \f$\mathcal{C}\f$ of the anisotropic stress power
 *     spectrum \f$k E_\Pi\f$.
 * 
 * This is the constant \f$\mathcal{C}\f$ defined in Eq. (46) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)",
 * \f[ 
 *     \mathcal{C} = \frac{16}{15} \int\limits_0^\infty d K \,
 *     \frac{\zeta_\mathrm{kin}^2(K)}{K^2} \,.
 * \f]
 * It describes the normalized amplitude of the anisotropic stress power
 * spectrum [cf. Eq. (44) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)"],
 * \f[
 *     k E_\Pi(\tau_1,\tau_2,k) = 2 \bar{w}^2 K^3 \left(
 *         \frac{E_\mathrm{kin}^*}{R_*}
 *     \right)^2 \mathcal{C} \zeta_\Pi(\tau_-, K)\,,
 * \f]
 * where \f$K = k R_*\f$.
 * 
 * @param kR Momentum \f$k\f$ times bubble separation \f$R_*\f$ at which the
 *     kinetic spectral shape is evaluated
 * @param zeta_kin Kinetic spectral shape \f$\zeta_\mathrm{kin}(k R_*)\f$
 * @return Normalized amplitude \f$\mathcal{C}\f$ of the anisotropic stress
 *     power spectrum \f$k E_\Pi\f$
 */
double anisotropic_stress_power_spectrum_amplitude(
    const std::vector<double>& kR,
    const std::vector<double>& zeta_kin
);

/**
 * @brief Calculate the nested momentum integral of the anisotropic stress UETC
 *     \f$E_\Pi\f$ or gravitational wave spectrum \f$\Omega_\mathrm{GW}\f$.
 * 
 * This function calculates the nested momentum (and time) integrals for
 * different choices of the time dependence function \f$\Delta\f$,
 * \f[
 *    \mathcal{I}(K, \Delta) =
 *     c_\bar{w} K^3 \int\limits_{-\infty}^{\infty} d P\,
 *     P^2 \bar{E}_\mathrm{kin}(P) \int\limits_{-1}^1 d z \, (1-z^2)^2
 *     \frac{\bar{E}_\mathrm{kin}(\tilde{P})}{\tilde{P}^4}
 *     \Delta(\delta \tau_\mathrm{fin}, K, P, \tilde{P})\,,
 * \f]
 * where \f$K = k R_*\f$, \f$P = p R_*\f$, \f$\tilde{P}^2 = K^2 + P^2
 * - 2 P K z\f$, and \f$\bar{E}_\mathrm{kin}(K) = R_*^{-1} E_\mathrm{kin}(k)\f$.
 * The choice of \f$\Delta\f$ is controlled by the \c universe and
 * \c freq_approx parameters of the @ref AnisotropicStressCorrelatorSettings
 * @p settings, and the prefactor \f$c_\bar{w}\f$ is set from `settings.c_w`. 
 * If `settings.universe` is set to `ExpansionHistory::NONE` (`settings.c_w` to
 * \f$2 \bar{w}^2\f$ where \f$\bar{w}\f$ is the mean enthalpy density), the
 * function returns the anisotropic stress (auto-correlation) power spectrum
 * \f$\mathcal{I}(K, 1) = k E_\Pi(k)\f$, where
 * \f$E_\Pi(k) = E_\Pi(\tau,\tau,k)\f$, i.e. Eq. (44) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
 * For `ExpansionHistory::RADIATION_DOMINATION` or
 * `ExpansionHistory::STATIC_UNIVERSE` (and `settings.c_w` = \f$3 \bar{w}^2\f$),
 * the function returns the gravitational wave spectrum \f$\mathcal{I}(K,\Delta)
 * = \Omega_\mathrm{GW}^*(k)\f$ at production, i.e. Eq. (51) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)" setting
 * \f$\mathrm{T}_\mathrm{GW}=1\f$, with \f$\Delta\f$ assuming an radiation
 * dominated or static universe, respectively, and using the approximation
 * specified in `settings.freq_approx`.
 * 
 * @param kR Momenta \f$k\f$ times bubble separation \f$R_*\f$ at which the
 *     anisotropic stress UETC is evaluated
 * @param E_kin_over_R Interpolating function for the scaled kinetic spectrum
 *     \f$R_*^{-1} E_\mathrm{kin}\f$ as a function of \f$\log(p R_*)\f$
 * @param settings @ref AnisotropicStressCorrelatorSettings object containing
 *     the settings for the momentum and time integration
 * @return Integrated anisotropic stress UETC \f$K^3 \bar{\zeta}(K, \Delta)\f$
 */
std::vector<double> anisotropic_stress_correlator_integration(
    const std::vector<double>& kR, 
    const utils::UniformLinearInterpolator& E_kin_over_R, 
    AnisotropicStressCorrelatorSettings settings
);

/**
 * @brief Anisotropic stress auto-correlation power spectrum \f$k E_\Pi\f$.
 *
 * This is the power spectrum of the anisotropic stress auto-correlation,
 * \f$k E_\Pi(k) = k E_\Pi(\tau,\tau,k)\f$, i.e. Eq. (47) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)" with
 * \f$\tau_-=0\f$.
 * \f[
 *     k E_\Pi(k) = 2 \bar{w}^2 K^3 \int\limits_{-\infty}^{\infty} d P\, 
 *     P^2 \bar{E}_\mathrm{kin}(P) \int\limits_{-1}^1 d z \, (1-z^2)^2
 *     \frac{\bar{E}_\mathrm{kin}(\tilde{P})}{\tilde{P}^4}\,,
 * \f]
 * where \f$K = k R_*\f$, i.e. the integration in
 * @ref anisotropic_stress_correlator_integration with \f$\Delta = 1\f$
 * (`settings.universe =`
 * `AnisotropicStressCorrelatorSettings::TimeDependence::None`)
 * and \f$c_\bar{w} = 2 \bar{w}^2\f$.
 * 
 * @param kR Momenta \f$k\f$ times bubble separation \f$R_*\f$ at which the
 *     spectrum is evaluated
 * @param E_kin_over_R Interpolating function for the scaled kinetic spectrum
 *     \f$R_*^{-1} E_\mathrm{kin}\f$ as a function of \f$\log(p R_*)\f$
 * @param w_bar Mean enthalpy density \f$\bar{w}\f$ (far ahead of the wall)
 *     (default: 1.0)
 * @param N_p Number of points in the integration over \f$p R_*\f$
 *     (default: 5000)
 * @param N_z Number of points in the integration over \f$z\f$ (default: 100) 
 * @return Anisotropic stress auto-correlator power spectrum
 *     \f$k E_\Pi(\tau,\tau,k)\f$
 */
std::vector<double> anisotropic_stress_power_spectrum(
    const std::vector<double>& kR,
    const utils::UniformLinearInterpolator& E_kin_over_R, 
    double w_bar=1.0,
    size_t N_p=5000,
    size_t N_z=100
);

/**
 * @brief Spectral shape of the anisotropic stress UETC power spectrum
 *     \f$(k R_*)^3 \zeta_\Pi(k R_*)\f$.
 * 
 * This is the combination
 * \f$K^3 \zeta_\Pi(K)\f$ with \f$\mathcal{C}\zeta_\Pi(K)\f$ defined in Eq. (47)
 * of  @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
 * It describes the spectral shape of the power spectrum of the anisotropic
 * stress UETC \f$k E_\Pi\f$,
 * \f[
 *     k E_\Pi(\tau_1,\tau_2,k) = 
 *     2 \bar{w}^2 \left(\frac{E_\mathrm{kin}^*}{R_*}\right)^2
 *     \mathcal{C} \times K^3 \zeta_\Pi(K)\,,
 * \f]
 * where \f$K = k R_*\f$ and the constant \f$\mathcal{C}\f$ can be calculated
 * with @ref anisotropic_stress_power_spectrum_amplitude
 * "anisotropic_stress_power_spectrum_amplitude(qR, zeta_kin)".
 * 
 * @param k_E_Pi Anisotropic stress auto-correlator \f$k E_\Pi\f$ 
 * @param script_C Amplitude parameter \f$\mathcal{C}\f$ of the anisotropic
 *     stress power spectrum
 * @param E_kin_max_over_R Maximal value \f$E_\mathrm{kin}^*\f$ of the kinetic
 *     spectrum \f$E_\mathrm{kin}\f$ normalized to the average bubble size
 *     \f$R_*\f$
 * @param w_bar Mean enthalpy density \f$\bar{w}\f$ (far ahead of the wall)
 *     (default: 1.0)
 * @return Spectral shape of the anisotropic stress UETC power spectrum
 *     \f$(k R_*)^3 \zeta_\Pi(k R_*)\f$
 */
std::vector<double> anisotropic_stress_power_spectral_shape(
    const std::vector<double>& k_E_Pi,
    double script_C,
    double E_kin_max_over_R,
    double w_bar=1.0
);

/**
 * @brief Approximated gravitational wave power spectrum
 *     \f$\Omega_\mathrm{GW}^*\f$ in the infinite duration limit.
 * 
 * This is the gravitational wave power spectrum \f$\Omega_\mathrm{GW}^*\f$ at
 * production [cf. @ref anisotropic_stress_correlator_integration] in the
 * inifinite source duration limit, Eq. (B3) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)" (with
 * \f$\mathcal{T}_\mathrm{GW}=1\f$). In the \f$k\gg 1/\delta\tau_\mathrm{fin}\f$
 * limit, the time integration factor \f$\Delta\f$ becomes a Dirac-\f$\delta\f$
 * distribution, simplifying the momentum integration to 
 * \f[
 *     \Omega_\mathrm{GW}^*(k) = c_\bar{w} \Upsilon(\tau_\mathrm{fin}) \,
 *     \mathcal{H}_* R_* \times \frac{\pi}{2 c_s}  
 *     \int\limits_\frac{(1-c_s) K}{2 c_s}^\frac{(1+c_s) K}{2 c_s} d P\, P 
 *     \bar{E}_\mathrm{kin}(P) (1-z^2)^2 
 *     \frac{\bar{E}_\mathrm{kin}(K/c_s-P)}{(K/c_s-P)^3}\,,
 * \f]
 * where \f$z=\frac{1}{c_s} - \frac{K(1-c_s^2)}{2 P c_s^2}\f$,
 * \f$\bar{E}_\mathrm{kin}(k R_*) = R_*^{-1} E_\mathrm{kin}(k)\f$, and we assume
 * \f$c_\bar{w} = 3 \Gamma^2\f$ with the mean adiabatic index in the stable
 * phase \f$\Gamma=\frac{\bar{w}}{\bar{e}}=\frac{4}{3}\f$ by default (although
 * this can be changed with @p c_w). The suppression factor \f$\Upsilon\f$ due
 * to the finite lifetime of the soundwave source is given by \f$
 *     \Upsilon =
 *     \frac{\mathcal{H}_* \delta\tau_\mathrm{fin}}
 *     {1 + \mathcal{H}_* \delta\tau_\mathrm{fin}}
 * \f$ in a radiation-dominated universe and
 * \f$\mathcal{H}_* \delta\tau_\mathrm{fin}\f$ in a static universe.
 * This corresponds to \f$ \mathcal{P}_\mathrm{GW}\f$ in Eq. (3.48) of
 * @ref ref_Hindmarsh_2019 "Hindmarsh and Hijazi (2019)".
 * 
 * @param kR Momenta \f$k\f$ times bubble separation \f$R_*\f$ at which the
 *     spectrum is evaluated
 * @param E_kin_over_R Interpolating function for the scaled kinetic spectrum
 *     \f$R_*^{-1} E_\mathrm{kin}\f$ as a function of \f$\log(p R_*)\f$
 * @param Upsilon_H_R Finite source lifetime suppression factor, \f$\Upsilon\f$,
 *     times the mean bubble separation normalized to the Hubble radius,
 *     \f$\mathcal{H}_* R_*\f$ 
 * @param N_p Number of points in the integration over \f$p R_*\f$
 *     (default: 5000)
 * @param c_w Prefactor \f$c_\bar{w}\f$ in
 *     @ref anisotropic_stress_correlator_integration, should be set to
 *     \f$c_\bar{w} = 3 \Gamma^2\f$ (default: 16.0/3.0)
 * @return GW power spectrum \f$\Omega_\mathrm{GW}^*\f$ at production in the
 *     \f$k \delta\tau_\mathrm{fin} \gg 1\f$ limit
 */
std::vector<double> infinite_duration_gw_power_spectrum(
    const std::vector<double>& kR,
    const utils::UniformLinearInterpolator& E_kin_over_R,
    double Upsilon_H_R,
    size_t N_p = 5000,
    double c_w = 16.0/3.0
);

/**
 * @brief Gravitational wave power spectrum \f$\Omega_\mathrm{GW}^*(k)\f$ at
 *     production.
 * 
 * This is Eq. (51) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)"
 * (with \f$\mathcal{T}_\mathrm{GW}=1\f$).
 * 
 * @param kR Momenta \f$k\f$ times bubble separation \f$R_*\f$ at which the
 *     spectrum is evaluated
 * @param E_kin_over_R Interpolating function for the scaled kinetic spectrum
 *     \f$R_*^{-1} E_\mathrm{kin}\f$ as a function of \f$\log(p R_*)\f$
 * @param H_delta_tau Source duration in units of the Hubble time,
 *     \f$\mathcal{H}_* \delta\tau_\mathrm{fin}\f$
 * @param H_R Mean bubble separation in units of the Hubble size,
 *     \f$\mathcal{H}_* R_*\f$
 * @param expansion Assumed @ref ExpansionHistory (default:
 *     \c ExpansionHistory::RADIATION_DOMINATION)
 * @param integration Method for calculating the momentum integrals in
 *     @ref anisotropic_stress_correlator_integration
 *     (default: \c MomentumIntegration::FULL)
 * @param Gamma Mean adiabatic index \f$\Gamma = \frac{\bar{w}}{\bar{e}}\f$ of
 *     the fluid in the stable phase (default: 4/3)
 * @param N_p Number of points in the integration over \f$p R_*\f$
 * @param N_z Number of points in the integration over \f$z\f$ 
 * @return GW power spectrum \f$\Omega_\mathrm{GW}^*(k)\f$ at production
 */
std::vector<double> gravitational_wave_power_spectrum(
    const std::vector<double>& kR,
    const utils::UniformLinearInterpolator& E_kin_over_R,
    double H_delta_tau, 
    double H_R,
    ExpansionHistory expansion=ExpansionHistory::RADIATION_DOMINATION,
    MomentumIntegration integration=MomentumIntegration::FULL,
    double Gamma=4.0/3.0,
    size_t N_p=5000,
    size_t N_z=100
);

/**
 * @brief Time integral function
 *     \f$\Delta(\delta \tau_\mathrm{fin},k,p,\tilde{p})\f$
 *     for a radiation dominated universe.
 * 
 * This is the function \f$\Delta\f$ that contains the integration over
 * conformal times \f$\tau_1\f$ and \f$\tau_2\f$ of the time dependence in the
 * Green's function and the stationary UETC [cf. Eq. (52) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)"], assuming a
 * radiation dominated universe.
 * \f[
 * \begin{aligned}
 *     \Delta(\delta \tau_\mathrm{fin}, k, p, \tilde{p}) &\equiv 
 *     \int\limits_{\tau_*}^{\tau_\mathrm{fin}} \frac{d \tau_1}{\tau_1} 
 *     \int\limits_{\tau_*}^{\tau_\mathrm{fin}} \frac{d \tau_2}{\tau_2} \, 
 *     \cos(p c_s \tau_-) \cos(\tilde{p} c_s \tau_-) \cos(k \tau_-) \\
 *     &= \sum\limits_{m,n=\pm 1} 
 *     \Delta_{mn}(\delta\tau_\mathrm{fin}, \hat{p}_{mn})\,,    
 * \end{aligned}
 * \f]
 * where \f$\tau_- = \tau_2-\tau_1\f$,
 * \f$\delta\tau_\mathrm{fin}=\tau_\mathrm{fin}-\tau_*\f$,
 * \f$\hat{p}_{mn} \equiv (p + m \tilde{p}) c_s + n k\f$,
 * and \f$\Delta_{mn}\f$ is defined in @ref Delta_mn.
 * 
 * We use the dimensionless variables \f$q \tau_* = \frac{q}{\mathcal{H}_*}\f$
 * for the wave numbers \f$q=k,p,\tilde{p}\f$ and \f$\mathcal{H}_*
 * \delta\tau_\mathrm{fin} = \frac{\delta\tau_\mathrm{fin}}{\tau_*}\f$ as
 * inputs. Note that, since we use \f$a_* = a(\tau_*) = 1\f$, we have
 * \f$H_* = \mathcal{H}_* = 1/\tau_*\f$, so that
 * \f$q \tau_* = \frac{q R_*}{H_* R_*}\f$.
 * 
 * @param H_delta_tau Duration of the gravitational wave source normalized to
 *     the Hubble time \f$\mathcal{H}_* \delta\tau_\mathrm{fin}\f$
 * @param k_tau Dimensionless wave number \f$k \tau_*\f$
 * @param p_tau Dimensionless wave number \f$p \tau_*\f$
 * @param pt_tau Dimensionless wave number \f$\tilde{p} \tau_*\f$
 * @return Value of the time integral \f$\Delta\f$
 */
double Delta_rad(double H_delta_tau, double k_tau, double p_tau, double pt_tau);

/**
 * @brief Intermediate function
 *     \f$\Delta_{mn}(\delta\tau_\mathrm{fin}, \hat{p}_{mn})\f$
 *     for the time integral in a radiation-dominated universe.
 * 
 * This is the intermediate function
 * \f$\Delta_{mn}(\delta\tau_\mathrm{fin}, \hat{p}_{mn})\f$ used in the
 * calculation of the time integral @ref Delta_rad [cf. Eq. (56) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)"],
 * \f[
 *     \Delta_{mn}(\delta\tau_\mathrm{fin}, \hat{p}_{mn}) = \frac{1}{4} \left\{
 *       \left[\mathrm{Ci}\left(\hat{p}_{mn} \tau_*
 *       (1+\mathcal{H}_*\delta\tau_\mathrm{fin})\right)
 *       - \mathrm{Ci}\left(\hat{p}_{mn} \tau_*\right)\right]^2
 *       + \left[\mathrm{Si}\left(\hat{p}_{mn} \tau_*
 *       (1+\mathcal{H}_*\delta\tau_\mathrm{fin})\right) 
 *       - \mathrm{Si}\left(\hat{p}_{mn} \tau_*\right)\right]^2
 *     \right\} \,.
 * \f]
 * We use the dimensionless variables
 * \f$\hat{p}_{mn} \tau_* = \frac{\hat{p}_{mn}}{\mathcal{H}_*}\f$ and
 * \f$\mathcal{H} \delta\tau_\mathrm{fin} =
 * \frac{\delta\tau_\mathrm{fin}}{\tau_*}\f$ as inputs.
 * 
 * @param H_delta_tau Duration of the gravitational wave source normalized to
 *     the Hubble time \f$\mathcal{H}_* \delta\tau_\mathrm{fin}\f$
 * @param pmn_tau Dimensionless wave number \f$\hat{p}_{mn} \tau_*\f$
 * @return Value of the intermediate time integral contribution
 *     \f$\Delta_{mn}\f$
 */
double Delta_mn(double H_delta_tau, double pmn_tau);

/**
 * @brief Time integral function \f$\Delta(\delta \tau_\mathrm{fin},k,p,
 *     \tilde{p})\f$ for a static universe.
 * 
 * This is the function \f$\Delta^\mathrm{static}\f$ that contains the
 * integration over conformal times \f$\tau_1\f$ and \f$\tau_2\f$ of the time
 * dependence in the Green's function and the stationary UETC [cf. Eq. (52) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)" with
 * \f$\frac{1}{\tau_1 \tau_2} \to \mathcal{H}_*^2\f$], assuming a static
 * universe.
 * \f[
 * \begin{aligned}
 *     \Delta(\delta \tau_\mathrm{fin}, k, p, \tilde{p}) &\equiv 
 *     \mathcal{H}_*^2 \int\limits_{\tau_*}^{\tau_\mathrm{fin}} d \tau_1
 * \int\limits_{\tau_*}^{\tau_\mathrm{fin}} d \tau_2 \, 
 *     \cos(p c_s \tau_-) \cos(\tilde{p} c_s \tau_-) \cos(k \tau_-) \\
 *     &= \sum\limits_{m,n=\pm 1}
 *     \Delta_{mn}^\mathrm{static}(\delta\tau_\mathrm{fin}, \hat{p}_{mn})\,,    
 * \end{aligned}
 * \f]
 * where \f$\tau_- = \tau_2-\tau_1\f$,
 * \f$\delta\tau_\mathrm{fin}=\tau_\mathrm{fin}-\tau_*\f$,
 * \f$\hat{p}_{mn} \equiv (p + m \tilde{p}) c_s + n k\f$,
 * and \f$\Delta_{mn}^\mathrm{static}\f$ is defined in @ref Delta_mn_static.
 * 
 * We use the dimensionless variables \f$q \tau_* = \frac{q}{\mathcal{H}_*}\f$
 * for the wave numbers \f$q=k,p,\tilde{p}\f$ and \f$\mathcal{H} \delta
 * \tau_\mathrm{fin} = \frac{\delta\tau_\mathrm{fin}}{\tau_*}\f$ as inputs.
 * Note that, since we use \f$a_* = a(\tau_*) = 1\f$, we have \f$H_* =
 * \mathcal{H}_* = 1/\tau_*\f$, so that \f$q \tau_* = \frac{q R_*}{H_* R_*}\f$.
 * 
 * @param H_delta_tau Duration of the gravitational wave source normalized to
 *     the Hubble time \f$\mathcal{H}_* \delta\tau_\mathrm{fin}\f$
 * @param k_tau Dimensionless wave number \f$k \tau_*\f$
 * @param p_tau Dimensionless wave number \f$p \tau_*\f$
 * @param pt_tau Dimensionless wave number \f$\tilde{p} \tau_*\f$
 * @return Value of the time integral \f$\Delta^\mathrm{flat}\f$ for a static
 *     universe
 */
double Delta_static(
    double H_delta_tau, double k_tau, double p_tau, double pt_tau
);

/**
 * @brief Intermediate function
 *     \f$\Delta_{mn}^\mathrm{flat}(\delta\tau_\mathrm{fin}, \hat{p}_{mn})\f$
 *     for the time integral in a static universe.
 * 
 * This is Eq. (59) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)",
 * \f[
 *     \Delta_{mn}^\mathrm{flat}(\delta\tau_\mathrm{fin}, \hat{p}_{mn}) =
 *     \frac{1-\cos(\hat{p}_{mn}\tau_*\,\mathcal{H}_* \delta\tau_\mathrm{fin})}
 *     {2 (\hat{p}_{mn}\tau_*)^2}\,.
 * \f]
 * 
 * We use the dimensionless variables
 * \f$\hat{p}_{mn} \tau_* = \frac{\hat{p}_{mn}}{\mathcal{H}_*}\f$ and
 * \f$\mathcal{H} \delta\tau_\mathrm{fin} =
 * \frac{\delta\tau_\mathrm{fin}}{\tau_*}\f$ as inputs.
 * 
 * @param H_delta_tau Duration of the gravitational wave source normalized to
 *     the Hubble time \f$\mathcal{H}_* \delta\tau_\mathrm{fin}\f$
 * @param pmn_tau Dimensionless wave number \f$\hat{p}_{mn} \tau_*\f$
 * @return Value of the time integral contribution
 *     \f$\Delta_{mn}^\mathrm{flat}\f$ neglecting the expansion of the Universe
 */
double Delta_mn_static(double H_delta_tau, double pmn_tau);

/**
 * @brief Time integral function in the low frequency limit
 *     \f$\Delta_0(\delta \tau_\mathrm{fin}, p,\tilde{p})
 *     = \lim_{k\to 0} \Delta(\delta \tau_\mathrm{fin},k,p,\tilde{p})\f$.
 * 
 * This is the low-frequency limit \f$\Delta_0\f$ of the integration over
 * conformal times \f$\tau_1\f$ and \f$\tau_2\f$ of the time dependence in the
 * Green's function and the stationary UETC [cf. Eq. (61) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)"],
 * assuming radiation domination.
 * \f[
 *     \Delta_0(\delta \tau_\mathrm{fin}, p) = \frac{1}{2} \left[
 *         \log^2(1+\mathcal{H}_* \delta\tau_\mathrm{fin})
 *       + \left[\mathrm{Ci}\left(
 *         2 p c_s \tau_* (1+\mathcal{H}_*\delta\tau_\mathrm{fin})
 *       \right) - \mathrm{Ci}\left(2 p c_s \tau_*\right)\right]^2
 *       + \left[\mathrm{Si}\left(
 *         2 p c_s \tau_* (1+\mathcal{H}_*\delta\tau_\mathrm{fin})
 *       \right) - \mathrm{Si}\left(2 p c_s \tau_*\right)\right]^2
 *     \right],
 * \f]
 * where \f$\delta\tau_\mathrm{fin}=\tau_\mathrm{fin}-\tau_*\f$.
 * 
 * We use the dimensionless variables \f$p \tau_* = \frac{p}{\mathcal{H}_*}\f$
 * and \f$\mathcal{H} \delta\tau_\mathrm{fin} =
 * \frac{\delta\tau_\mathrm{fin}}{\tau_*}\f$ as inputs. Note that, since we use
 * \f$a_* = a(\tau_*) = 1\f$, we have \f$H_* = \mathcal{H}_* = 1/\tau_*\f$, so
 * that \f$p \tau_* = \frac{p R_*}{H_* R_*}\f$.
 * 
 * @param H_delta_tau Duration of the gravitational wave source normalized to
 *     the Hubble time \f$\mathcal{H}_* \delta\tau_\mathrm{fin}\f$
 * @param p_tau Dimensionless wave number \f$p \tau_*\f$
 * @return Value of the time integral \f$\Delta_0\f$ in the low frequency limit
 */
double Delta_low_freq_rad(double H_delta_tau, double p_tau);

/**
 * @brief Time integral function in the low frequency limit in a static
 *     universe,
 *     \f$\Delta_0^\mathrm{flat}(\delta \tau_\mathrm{fin}, p) = \lim_{k\to 0}
 *     \Delta^\mathrm{flat}(\delta \tau_\mathrm{fin},k,p,\tilde{p})\f$.
 * 
 * This is the low-frequency limit \f$\Delta_0^\mathrm{flat}\f$ of the
 * integration over conformal times \f$\tau_1\f$ and \f$\tau_2\f$ of the time
 * dependence in the Green's function and the stationary UETC [cf. Eq. (61) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)"], assuming a
 * static universe.
 * \f[
 *     \Delta_0^\mathrm{flat}(\delta \tau_\mathrm{fin}, p)
 *     = \frac{1}{2} \left(\mathcal{H}_* \delta\tau_\mathrm{fin}\right)^2 \left[
 *         1 + \mathrm{sinc}^2\left(p c_s \delta\tau_\mathrm{fin}\right)
 *     \right],
 * \f]
 * where \f$\delta\tau_\mathrm{fin}=\tau_\mathrm{fin}-\tau_*\f$.
 * 
 * We use the dimensionless variables \f$p \tau_* = \frac{p}{\mathcal{H}_*}\f$
 * and \f$\mathcal{H} \delta\tau_\mathrm{fin} =
 * \frac{\delta\tau_\mathrm{fin}}{\tau_*}\f$ as inputs. Note that, since we use
 * \f$a_* = a(\tau_*) = 1\f$, we have \f$H_* = \mathcal{H}_* = 1/\tau_*\f$, so
 * that \f$q \tau_* = \frac{q R_*}{H_* R_*}\f$.
 * 
 * @param H_delta_tau Duration of the gravitational wave source normalized to
 *     the Hubble time \f$\mathcal{H}_* \delta\tau_\mathrm{fin}\f$
 * @param p_tau Dimensionless wave number \f$p \tau_*\f$
 * @return Value of the time integral \f$\Delta_0^\mathrm{flat}\f$ in the low
 *     frequency limit in a static universe
 */
double Delta_low_freq_static(double H_delta_tau, double p_tau);

/**
 * @brief Weighted average \f$\tilde{\Delta}_0\f$ of \f$\Delta_0\f$ weighted
 *     with the kinetic spectrum \f$E_\mathrm{kin}\f$.
 * 
 * This is the weighted average of the low-frequency integral function
 * \f$\Delta_0\f$ with the kinetic spectrum \f$E_\mathrm{kin}\f$,
 * \f[
 *     \tilde{\Delta}_0(\delta\tau_\mathrm{fin}, R_*) = \frac{
 *         \int\limits_0^\infty d(k R_*) \,\frac{E_\mathrm{kin}^2(k)}{(k R_*)^2}
 *          \Delta_0(\delta\tau_\mathrm{fin}, k)
 *     }{
 *         \int\limits_0^\infty d(k R_*) \,\frac{E_\mathrm{kin}^2(k)}{(k R_*)^2}
 *     }\,,
 * \f]
 * parametrizing the dependence of the GW amplitude on the source duration
 * \f$\delta\tau_\mathrm{fin}\f$ [cf. Eq. (63) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)"].
 * 
 * @param kR Momenta \f$k\f$ times bubble separation \f$R_*\f$
 * @param E_kin_over_R Scaled kinetic spectrum
 *     \f$\bar{E}_\mathrm{kin}(k R_*) = R_*^{-1} E_\mathrm{kin}(k)\f$
 * @param H_delta_tau Duration of the gravitational wave source normalized to
 *     the Hubble time \f$\mathcal{H}_* \delta\tau_\mathrm{fin}\f$
 * @param H_R Mean bubble separation \f$\mathcal{H}_* R_*\f$ normalized to the
 *     Hubble radius
 * @param expansion Assumed @ref ExpansionHistory (default:
 *     ExpansionHistory::RADIATION_DOMINATION)
 * @return Weighted average source lifetime dependence
 *     \f$\tilde{\Delta}_0(\delta\tau_\mathrm{fin}, R_*)\f$
 */
double weighted_average_Delta_0(
    const std::vector<double>& kR,
    const std::vector<double>& E_kin_over_R,
    double H_delta_tau,
    double H_R, 
    ExpansionHistory expansion=ExpansionHistory::RADIATION_DOMINATION
);

/**
 * @brief Normalized gravitational wave spectrum
 *     \f$(k R_*)^3 \zeta_\mathrm{GW}(k R_*)\f$.
 * 
 * This is the normalized gravitational wave spectrum
 * \f[
 *     K^3 \zeta_\mathrm{GW}(K) = \left(\frac{R_*}{E_\mathrm{kin}^*}\right)^2
 *     \frac{
 *       \Omega_\mathrm{GW}(\delta\tau_\mathrm{fin}, R_*, K)
 *     }{
 *       3 \Gamma^2 \mathcal{T}_\mathrm{GW} \mathcal{C}
 *       \tilde{\Delta}_0(\delta\tau_\mathrm{fin}, R_*)
 *     }
 * \f]
 * in Eq. (93) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
 * 
 * @param Omega_GW Gravitational wave power spectrum
 *     \f$\Omega_\mathrm{GW}^*(k)\f$ at production,
 *     cf. @ref gravitational_wave_power_spectrum
 * @param script_C Amplitude parameter \f$\mathcal{C}\f$ of the anisotropic
 *     stress power spectrum,
 *     cf. @ref anisotropic_stress_power_spectrum_amplitude
 * @param E_kin_max_over_R Maximal value \f$E_\mathrm{kin}^*\f$ of the kinetic
 *     spectrum \f$E_\mathrm{kin}\f$ normalized to the average bubble size
 *     \f$R_*\f$
 * @param Delta_tilde_0 Weighted average source lifetime dependence
 *     \f$\tilde{\Delta}_0\f$, cf. @ref weighted_average_Delta_0
 * @param Gamma Mean adiabatic index \f$\Gamma = \frac{w}{\bar{e}}\f$ of the
 *     fluid in the stable phase (default: 4/3)
 * @return Normalized gravitational wave spectrum \f$K^3 \zeta_\mathrm{GW}\f$
 */
std::vector<double> normalized_GW_spectrum(
    const std::vector<double>& Omega_GW, 
    double script_C,
    double E_kin_max_over_R,
    double Delta_tilde_0,
    double Gamma = 4.0/3.0
);

/**
 * @brief Spectral modification function
 *     \f$\tilde{\Delta} = \zeta_\mathrm{GW}/zeta_{\Pi}\f$.
 * 
 * This functions parametrized the modification of the spectral shape of the
 * graviational wave spectrum due to the finite source duration, Eq. (95) of
 * @ref ref_Roper_Pol_2024 "Roper Pol, Procacci and Caprini (2024)".
 * 
 * @param zeta_GW Normalized gravitational wave spectrum \f$\zeta_\mathrm{GW}\f$
 *     (or \f$K^3 \zeta_\mathrm{GW}\f$) [cf. @ref normalized_GW_spectrum]
 * @param zeta_Pi Normalized anisotropic stress spectrum \f$\zeta_\Pi\f$ (or
 *     \f$K^3 \zeta_\Pi\f$) [cf. @ref anisotropic_stress_power_spectral_shape]
 * @return Spectral modification function
 *     \f$\tilde{\Delta} = \zeta_\mathrm{GW}/zeta_{\Pi}\f$
 */
std::vector<double> spectral_modification_function(
    const std::vector<double>& zeta_GW, const std::vector<double>& zeta_Pi
);

}; // sound_shell_model