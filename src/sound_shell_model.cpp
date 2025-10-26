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

#include "inverse_pt/sound_shell_model.hpp"

#include "inverse_pt/constants.hpp"
#include "inverse_pt/settings.hpp"
#include "inverse_pt/utils.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <omp.h>

#include <gsl/gsl_sf_expint.h>

using vecd = std::vector<double>;
using namespace inverse_pt::utils;
using namespace inverse_pt::constants;
using namespace inverse_pt::settings;

__attribute__((always_inline))
inline double Ci(double x) {
    return gsl_sf_Ci(x);
}

__attribute__((always_inline))
inline double Si(double x) {
    return gsl_sf_Si(x);
}

namespace inverse_pt::sound_shell_model {

std::tuple<vecd, vecd, vecd> calculate_sine_transforms(
    const vecd& chi_vec, 
    const vecd& xi, 
    const vecd& v, 
    const vecd& w, 
    const vecd& e, 
    size_t idx_wall, 
    size_t idx_sh
) {
    const size_t N = xi.size();
    if ( N != v.size() || N != w.size() || N != e.size() || N < 2 ) {
        throw std::invalid_argument(
            "Fluid profile vectors must have same size (>2). Got sizes " 
            + std::to_string(N) + ", " + std::to_string(v.size()) + "," 
            + std::to_string(w.size()) + " and " + std::to_string(e.size())
            + " for xi, v, w and e, respectively."
        );
    }
    // mean energy density (far in front of the bubble)
    const double ebar = e[N-1]; 
    // mean enthalpy density (far in front of the bubble)
    const double w_N = w[N-1]; 

    const size_t N_chi = chi_vec.size();
    vecd A2_vec(N_chi), df_vec(N_chi), l_vec(N_chi);
    // loop over z-values
    #pragma omp parallel for // only outer loop parallelized!
    for ( size_t j = 0; j < N_chi; ++j ) {
        double chi = chi_vec[j];
        double prefactor = 4.0 * Pi / chi / chi / chi;

        // precompute sine and cosine
        std::vector<double> x(N), sin_x(N), cos_x(N);
        #pragma omp simd
        for ( size_t i = 0; i < N; ++i ) {
            x[i] = chi * xi[i];
            sin_x[i] = std::sin(x[i]);
            cos_x[i] = std::cos(x[i]);
        }

        // contribution of integral boundaries to df/dz and l
        auto df_jump = [x, v, cos_x, sin_x](size_t idx) {
            return v[idx] * (2.0 * cos_x[idx] + x[idx] * sin_x[idx]); 
        };
        auto l_jump = [x, e, cos_x, sin_x, ebar, w_N](size_t idx) {
            return (e[idx] - ebar) / w_N * (sin_x[idx] - x[idx] * cos_x[idx]);
        };

        // compute df, l and A-squared
        double df = (
            df_jump(N-1) - df_jump(idx_wall) + df_jump(idx_wall-1) - df_jump(0)
        );
        double l = - l_jump(idx_wall) + l_jump(idx_wall-1) - l_jump(0);
        if ( idx_sh != NO_INDEX ) {
            df += df_jump(idx_sh-1) - df_jump(idx_sh);
            l += l_jump(idx_sh-1) - l_jump(idx_sh);
        }
        #pragma omp simd reduction(+:df, l)
        for ( size_t i = 0; i < N - 1; ++i ) {
            const double over_dx = (
                std::abs(xi[i+1] - xi[i]) >= jump_tolerance ?
                1.0 / (x[i+1] - x[i]) : 0.0
            );
            const double slope_df = (v[i+1] - v[i]) * over_dx;
            const double slope_integral_df = (
                x[i+1] * cos_x[i+1] - x[i] * cos_x[i] 
                + 3.0 * (sin_x[i] - sin_x[i+1])
            );
            const double slope_l = (e[i+1] - e[i]) * over_dx / w_N;
            const double slope_integral_l = (
                x[i+1] * sin_x[i+1] - x[i] * sin_x[i] 
                + 2.0 * cos_x[i+1] - 2.0 * cos_x[i]
            );
            df += slope_df * slope_integral_df;
            l += slope_l * slope_integral_l;
        }
        df *= prefactor;
        l *= prefactor;   
        double A2 = 0.25 * (df * df + cs2 * l * l);

        // save results
        A2_vec[j] = A2;
        df_vec[j] = df;
        l_vec[j] = l; 
    }

    return {A2_vec, df_vec, l_vec};
}

vecd scaled_kinetic_spectrum(
    const vecd& kR, 
    const vecd& chi_vec, 
    const vecd& A2_vec, 
    double xi_wall, 
    NucleationHistory nuc_hist
) {
    size_t N_k = kR.size(), N_chi = chi_vec.size();
    if ( N_chi != A2_vec.size() || N_chi < 1 ) { throw std::invalid_argument(
        "Vectors chi and A2 must have the same size and cannot be empty."
        " Got sizes " + std::to_string(N_chi) + " and "
        + std::to_string(A2_vec.size()) + "."
    ); }

    // third root of pi
    constexpr double cbrt_pi = 1.464591887561523; 
    // beta * R = (8 pi)^(1/3) v_w
    const double beta_R = 2.0 * cbrt_pi * xi_wall; 
    // 1/(beta*R)^6/(2 pi^2) = 1/(128 pi^4)/v_w^6
    const double prefactor = (
        0.0078125 / std::pow(inverse_pt::constants::Pi, 4) 
                  / std::pow(xi_wall, 6)
    ); 

    std::vector<double> E_kin_over_R(N_k);
    #pragma omp parallel for
    for( size_t i = 0; i < N_k; ++i ) {

        // precompute integrand
        std::vector<double> integrand(N_chi);
        if ( nuc_hist == NucleationHistory::EXPONENTIAL ) {
            #pragma omp simd
            for ( size_t j = 0; j < N_chi; ++j )
            {
                const double T = beta_R * chi_vec[j] / kR[i];
                integrand[j] = T*T*T*T*T*T * nu_exp(T) * A2_vec[j];
            }
        } else if ( nuc_hist == NucleationHistory::SIMULTANEOUS ) {
            #pragma omp simd
            for ( size_t j = 0; j < N_chi; ++j )
            {
                const double T = beta_R * chi_vec[j] / kR[i];
                integrand[j] = T*T*T*T*T*T * nu_sim(T) * A2_vec[j];
            }
        } else {
            throw std::runtime_error("Cannot identify nucleation history.");
        }

        // calculate integral
        E_kin_over_R[i] = (
            prefactor * beta_R * kR[i] 
            * utils::integrate_trapez(chi_vec, integrand)
        ); // kR^2 * beta_R/kR to change chi integration to T integration
    }

    return E_kin_over_R;
}

vecd kinetic_spectral_shape(const vecd& E_kin_over_R) {
    double E_star_over_R = *std::max_element(
        E_kin_over_R.begin(), E_kin_over_R.end()
    );
    return utils::rescale_vector(E_kin_over_R, 1.0/E_star_over_R);
}

double mean_squared_fluid_velocity(const vecd& kR, const vecd& E_kin_over_R) { 
    return 2.0 * utils::integrate_log_trapez(kR, E_kin_over_R); 
}

double mean_squared_fluid_velocity(
    const vecd& chi_vec, const vecd& A2_vec, double xi_wall
) { 
    size_t N_chi = chi_vec.size();
    if ( N_chi != A2_vec.size() || N_chi < 1 ) { throw std::invalid_argument(
        "Vectors chi_vec and A2_vec must have the same size and cannot be"
        " empty. Got sizes " + std::to_string(N_chi) + " and " 
        + std::to_string(A2_vec.size()) + "."
    ); }
    const double prefactor = 3.0 / std::pow(2.0 * Pi * xi_wall, 3);
    std::vector<double> integrand(N_chi);
    #pragma omp simd 
    for ( size_t i = 0; i < N_chi; ++i ) { 
        integrand[i] = chi_vec[i]*chi_vec[i]*A2_vec[i];
    }
    return prefactor * utils::integrate_log_trapez(chi_vec, integrand);
}

double kinetic_energy_density_parameter(double Uf2, double Gamma) { 
    return Gamma * Uf2; 
}

double kinetic_energy_density_parameter(
    const vecd& kR, const vecd& E_kin_over_R, double Gamma
) {
    return 2.0 * Gamma * utils::integrate_log_trapez(kR, E_kin_over_R);
}

double kinetic_energy_density_parameter(
    double xi_wall, const vecd& chi_vec, const vecd& A2_vec, double Gamma
) {
    return Gamma * mean_squared_fluid_velocity(chi_vec, A2_vec, xi_wall);
}

double anisotropic_stress_power_spectrum_amplitude(
    const vecd& kR, const vecd& zeta_kin
) {
    size_t N_k = kR.size();
    if ( N_k != zeta_kin.size() || N_k < 1 ) { throw std::invalid_argument(
        "Vectors kR and zeta_kin must have the same size and cannot be empty."
        " Got sizes " + std::to_string(N_k) + " and " 
        + std::to_string(zeta_kin.size()) + "."
    ); }
    std::vector<double> integrand(N_k);
    #pragma omp simd 
    for ( size_t i = 0; i < N_k; ++i ) { 
        integrand[i] = zeta_kin[i]*zeta_kin[i]/kR[i]/kR[i];
    }
    return 16.0/15.0 * utils::integrate_log_trapez(kR, integrand);
}

vecd anisotropic_stress_correlator_integration(
    const vecd& kR,
    const utils::UniformLinearInterpolator& E_kin_over_R,
    AnisotropicStressCorrelatorSettings settings
) { 
    if ( settings.universe == ExpansionHistory::NONE ) {
        if ( 
            settings.freq_approx 
            != AnisotropicStressCorrelatorSettings::MomentumIntegration::FULL 
        ) { throw std::invalid_argument(
            "Low-frequency and infinite-duration appoximations cannot be used"
            " for the autocorrelation E_Pi without time integration.\n"
            "Either set 'universe' to 'RADIATION_DOMINATION' or"
            " 'STATIC_UNIVERSE' or use 'freq_approx = Full'."
        ); }
    } else { if ( settings.H_delta_tau <= 0.0 || settings.H_R <= 0.0 ) { 
        throw std::invalid_argument(
            "H_delta_tau and H_R must be positive for time integration. Got " 
            + std::to_string(settings.H_delta_tau) + " and " 
            + std::to_string(settings.H_R) + "."
        ); 
    } }
    if ( 
        settings.freq_approx == AnisotropicStressCorrelatorSettings::
            MomentumIntegration::INFINITE_DURATION 
    ) {
        double Upsilon_H_R;
        if ( settings.universe == ExpansionHistory::RADIATION_DOMINATION ) {
            Upsilon_H_R = settings.H_R * settings.H_delta_tau/(
                1.0 + settings.H_delta_tau
            );
        } else if ( settings.universe == ExpansionHistory::STATIC_UNIVERSE ) {
            Upsilon_H_R = settings.H_R * settings.H_delta_tau;
        } else { throw std::invalid_argument("Invalid time dependence type."); }
        return infinite_duration_gw_power_spectrum(
            kR, E_kin_over_R, Upsilon_H_R, settings.N_p, settings.c_w
        );
    }
    size_t N_k = kR.size();
    size_t N_z = settings.N_z;
    size_t N_p = settings.N_p;
    double H_delta_tau = settings.H_delta_tau;
    double H_R = settings.H_R;
    double c_w = settings.c_w;
    std::vector<double> z = lin_space(-1.0, 1.0, N_z);
    std::vector<double> pR = log_space(
        std::exp(E_kin_over_R.get_xmin()), 
        std::exp(E_kin_over_R.get_xmax()), 
        N_p
    );
    std::vector<double> log_pR = lin_space(
        E_kin_over_R.get_xmin(), E_kin_over_R.get_xmax(), N_p
    );
    double dz = 2.0 / (N_z - 1.0);
    double d_log_pR = (
        E_kin_over_R.get_xmax() - E_kin_over_R.get_xmin()
    ) / (N_p - 1.0);
    // vector for storing the result
    std::vector<double> result(N_k, 0.0); 
    // main loop over k
    #pragma omp parallel shared(E_kin_over_R)
    {
        std::vector<double> ptR2(N_p);
        std::vector<double> Delta_vals(N_p, 1.0);
        // for low-frequency limit, precompute Delta here
        if ( 
            settings.freq_approx == AnisotropicStressCorrelatorSettings::
                MomentumIntegration::LOW_FREQUENCY
        ) {
            if ( settings.universe == ExpansionHistory::RADIATION_DOMINATION ) {
                //#pragma omp simd
                for ( size_t k = 0; k < N_p; ++k ) { 
                    Delta_vals[k] = Delta_low_freq_rad(H_delta_tau, pR[k]/H_R); 
                }
            } else if ( settings.universe==ExpansionHistory::STATIC_UNIVERSE ) {
                #pragma omp simd
                for ( size_t k = 0; k < N_p; ++k ) { 
                    Delta_vals[k] = Delta_low_freq_static(
                        H_delta_tau, pR[k]/H_R
                    );
                }
            }
        }
        #pragma omp for
        for ( size_t i = 0; i < N_k; ++i ) {
            double integral = 0.0;
            for ( size_t j = 0; j < N_z; ++j ) {
                double z_factor = (
                    (1.0 - z[j]*z[j])*(1.0 - z[j]*z[j]) 
                    * ((j==0||j==N_z-1) ? 0.5 : 1.0)
                );
                // precompute p-tilde
                #pragma omp simd
                for ( size_t k = 0; k < N_p; ++k ) {
                    ptR2[k] = kR[i] * kR[i] + pR[k] * pR[k] 
                              - 2.0 * kR[i] * pR[k] * z[j];
                }
                // For 'FULL' momentum integration, need to precompute Delta
                // here
                if ( 
                    settings.freq_approx == 
                    AnisotropicStressCorrelatorSettings::
                    MomentumIntegration::FULL 
                ) {
                    if ( 
                        settings.universe == 
                        ExpansionHistory::RADIATION_DOMINATION 
                    ) {
                        for ( size_t k = 0; k < N_p; ++k ) {
                            Delta_vals[k] = Delta_rad(
                                H_delta_tau, 
                                kR[i]/H_R, 
                                pR[k]/H_R, 
                                std::sqrt(ptR2[k])/H_R
                            );
                        }
                    } else if ( 
                        settings.universe == ExpansionHistory::STATIC_UNIVERSE 
                    ) {
                        #pragma omp simd
                        for ( size_t k = 0; k < N_p; ++k ) {
                            Delta_vals[k] = Delta_static(
                                H_delta_tau, 
                                kR[i]/H_R, 
                                pR[k]/H_R, 
                                std::sqrt(ptR2[k])/H_R
                            );
                        }
                    }
                }
                // perform integration
                double inner_integral = 0.5 * (
                    E_kin_over_R(log_pR[0]) * pR[0]*pR[0]*pR[0] 
                    * E_kin_over_R(0.5*std::log(ptR2[0]))/ptR2[0]/ptR2[0] 
                    * Delta_vals[0]
                    + E_kin_over_R(log_pR[N_p-1]) 
                    * pR[N_p-1]*pR[N_p-1]*pR[N_p-1] 
                    * E_kin_over_R(0.5*std::log(ptR2[N_p-1]))
                    /ptR2[N_p-1]/ptR2[N_p-1] * Delta_vals[N_p-1]
                );
                #pragma omp simd reduction(+:inner_integral)
                for ( size_t k = 1; k < N_p-1; ++k ) {
                    inner_integral +=
                        E_kin_over_R(log_pR[k]) * pR[k]*pR[k]*pR[k] 
                        * E_kin_over_R(0.5*std::log(ptR2[k]))/ptR2[k]/ptR2[k] 
                        * Delta_vals[k];
                }
                integral += inner_integral * z_factor;
            }
            result[i] += c_w * kR[i]*kR[i]*kR[i] * integral * d_log_pR * dz;
        }
    }
    return result;
};

vecd anisotropic_stress_power_spectrum(
    const vecd& kR, 
    const utils::UniformLinearInterpolator& E_kin_over_R, 
    double w_bar, 
    size_t N_p, 
    size_t N_z
) {
    return anisotropic_stress_correlator_integration(
        kR, 
        E_kin_over_R, 
        AnisotropicStressCorrelatorSettings(N_p, N_z, 2.0*w_bar*w_bar)
    );
}

vecd anisotropic_stress_power_spectral_shape(
    const vecd& k_E_Pi, double script_C, double E_kin_max_over_R, double w_bar
) {
    return utils::rescale_vector(
        k_E_Pi, 0.5/script_C/E_kin_max_over_R/E_kin_max_over_R/w_bar/w_bar
    );
}

vecd infinite_duration_gw_power_spectrum(
    const vecd& kR,
    const utils::UniformLinearInterpolator& E_kin_over_R,
    double Upsilon_H_R,
    size_t N_p,
    double c_w
) {
    size_t N_k = kR.size();
    std::vector<double> result(N_k); 
    double log_p_over_k_min = std::log(0.5 * (1.0 - cs)/cs);
    double log_p_over_k_max = std::log(0.5 * (1.0 + cs)/cs);
    #pragma omp parallel for
    for ( size_t i = 0; i < N_k; ++i ) {
        double integral = 0.0;
        double log_pR_min = std::max(
            E_kin_over_R.get_xmin(), std::log(kR[i]) + log_p_over_k_min
        );
        double log_pR_max = std::min(
            E_kin_over_R.get_xmax(), std::log(kR[i]) + log_p_over_k_max
        );
        double d_log_pR = (log_pR_max - log_pR_min) / (N_p - 1.0);
        if ( d_log_pR <= 0.0 ) { result[i] = 0.0; continue; }
        #pragma omp simd reduction(+:integral)
        for ( size_t j = 0; j < N_p; ++j ) {
            double C = (j==0||j==N_p-1) ? 0.5 : 1.0;
            double log_pR = log_pR_min + j*d_log_pR;
            double pR = std::exp(log_pR);
            double ptR = kR[i]/cs - pR;
            double z = 1.0/cs - 0.5*kR[i]*(1.0-cs2)/pR/cs2;
            double integrand = (
                pR*pR/ptR/ptR/ptR * (1.0-z*z)*(1.0-z*z) 
                * E_kin_over_R(log_pR) * E_kin_over_R(std::log(ptR))
            );
            integral += C * integrand * d_log_pR; 
        }
        result[i] = c_w * 0.5*Pi * Upsilon_H_R * kR[i]*kR[i] / cs * integral;
    }
    return result;
}

vecd gravitational_wave_power_spectrum(
    const vecd& kR,
    const utils::UniformLinearInterpolator& E_kin_over_R,
    double H_delta_tau, double H_R,
    ExpansionHistory expansion,
    MomentumIntegration integration,
    double Gamma, size_t N_p, size_t N_z
) {
    if ( expansion == ExpansionHistory::NONE ) { throw std::invalid_argument(
        "ExpansionHistory::NONE is not a valid expansion history for"
        " calculating the GW spectrum. Please use ExpansionHistory::"
        "RadiationDomination or ExpansionHistory::StaticUniverse instead."
    );}
    return anisotropic_stress_correlator_integration(
        kR,
        E_kin_over_R,
        AnisotropicStressCorrelatorSettings(
            N_p, N_z, H_delta_tau, H_R, expansion, integration, 3.0*Gamma*Gamma
        )
    );
}

double Delta_mn(double H_delta_tau, double pmn_tau) {
    double DeltaCi = Ci(pmn_tau * (1.0 + H_delta_tau)) - Ci(pmn_tau);
    double DeltaSi = Si(pmn_tau * (1.0 + H_delta_tau)) - Si(pmn_tau);
    double result = 0.25 * (DeltaCi * DeltaCi + DeltaSi * DeltaSi);
    return result;
}

double Delta_rad(
    double H_delta_tau, double k_tau, double p_tau, double pt_tau
) {
    return Delta_mn(H_delta_tau, std::abs((p_tau - pt_tau) * cs - k_tau)) +
           Delta_mn(H_delta_tau, std::abs((p_tau - pt_tau) * cs + k_tau)) +
           Delta_mn(H_delta_tau, std::abs((p_tau + pt_tau) * cs - k_tau)) +
           Delta_mn(H_delta_tau, std::abs((p_tau + pt_tau) * cs + k_tau));
}

double Delta_static(
    double H_delta_tau, double k_tau, double p_tau, double pt_tau
) {
    return Delta_mn_static(H_delta_tau, (p_tau - pt_tau) * cs - k_tau) +
           Delta_mn_static(H_delta_tau, (p_tau - pt_tau) * cs + k_tau) +
           Delta_mn_static(H_delta_tau, (p_tau + pt_tau) * cs - k_tau) +
           Delta_mn_static(H_delta_tau, (p_tau + pt_tau) * cs + k_tau);
}

double Delta_mn_static(double H_delta_tau, double pmn_tau) {
    return 0.5 * (1.0 - std::cos(pmn_tau*H_delta_tau)) /pmn_tau/pmn_tau;
}

double Delta_low_freq_rad(double H_delta_tau, double p_tau) {
    double log_H_tau = std::log(1.0 + H_delta_tau);
    double pmn_tau = 2.*p_tau*cs;
    double DeltaCi = Ci(pmn_tau*(1.0+H_delta_tau)) - Ci(pmn_tau);
    double DeltaSi = Si(pmn_tau*(1.0+H_delta_tau)) - Si(pmn_tau);
    return 0.5 * (log_H_tau*log_H_tau + DeltaCi*DeltaCi + DeltaSi*DeltaSi);
}

double Delta_low_freq_static(double H_delta_tau, double p_tau) {
    double x = p_tau*cs*H_delta_tau;
    double sinc = (x==0.0) ? 1.0 : std::sin(x) / x;
    return 0.5 * H_delta_tau*H_delta_tau * (1.0 + sinc*sinc);
}

double weighted_average_Delta_0(
    const vecd& kR, 
    const vecd& E_kin_over_R, 
    double H_delta_tau, 
    double H_R, 
    ExpansionHistory expansion
) {
    if ( expansion == ExpansionHistory::NONE ) { throw std::invalid_argument(
        "ExpansionHistory::NONE is not a valid expansion history for"
        " calculating the source duration dependence. Please use"
        " ExpansionHistory::RadiationDomination or ExpansionHistory::"
        "StaticUniverse instead."
    );}
    size_t N_k = kR.size();
    if ( N_k != E_kin_over_R.size() || N_k < 1 ) { throw std::invalid_argument(
        "Vectors kR and E_kin_over_R must have the same size and cannot be"
        " empty. Got sizes " + std::to_string(N_k) + " and " 
        + std::to_string(E_kin_over_R.size()) + "."
    ); }
    std::vector<double> log_kR(N_k), integrand(N_k);
    #pragma omp simd 
    for ( size_t i = 0; i < N_k; ++i ) {
        log_kR[i] = std::log(kR[i]);
        integrand[i] = E_kin_over_R[i]*E_kin_over_R[i]/kR[i];
    }
    double denominator = utils::integrate_trapez(integrand, log_kR);
    if ( expansion == ExpansionHistory::RADIATION_DOMINATION ) {
        #pragma omp parallel for
        for ( size_t i = 0; i < N_k; ++i ) { 
            integrand[i] *= Delta_low_freq_rad(H_delta_tau, kR[i]/H_R); 
        }
    } else {
        #pragma omp simd
        for ( size_t i = 0; i < N_k; ++i ) { 
            integrand[i] *= Delta_low_freq_static(H_delta_tau, kR[i]/H_R); 
        }
    }
    return utils::integrate_trapez(integrand, log_kR) / denominator;
}

vecd normalized_GW_spectrum(
    const vecd& Omega_GW, 
    double script_C, 
    double E_kin_max_over_R, 
    double Delta_tilde_0, 
    double Gamma
) {
    return utils::rescale_vector(
        Omega_GW, 
        1.0 / ( 
            3.0 * Gamma*Gamma * script_C 
            * E_kin_max_over_R * E_kin_max_over_R 
            * Delta_tilde_0
        )
    );
}

vecd spectral_modification_function(const vecd& zeta_GW, const vecd& zeta_Pi) {
    size_t N = zeta_GW.size();
    if ( N != zeta_Pi.size() || N < 1 ) { throw std::invalid_argument(
        "Vectors zeta_GW and zeta_Pi must have the same size and cannot be"
        " empty. Got sizes " + std::to_string(N) + " and " 
        + std::to_string(zeta_Pi.size()) + "."
    ); } 
    std::vector<double> result(N);
    #pragma omp simd
    for ( size_t i = 0; i < N; ++i ) { result[i] = zeta_GW[i]/zeta_Pi[i]; }
    return result;
}

} // namespace sound_shell_model
