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

#include "inverse_pt/sound_shell_spectrum.hpp"

#include "inverse_pt/fluid_profile.hpp"
#include "inverse_pt/sound_shell_model.hpp"
#include "inverse_pt/utils.hpp"
#include "inverse_pt/constants.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include <omp.h>

using NucleationHistory = inverse_pt::sound_shell_model::NucleationHistory;
using ExpansionHistory = inverse_pt::sound_shell_model::ExpansionHistory;

namespace inverse_pt {

SoundShellSpectrum::SoundShellSpectrum(
    const FluidProfile& profile, 
    double R_H,
    NucleationHistory nhist,
    ExpansionHistory exhist
) 
: fluid_profile(profile), nuc_hist(nhist), exp_hist(exhist), H_R(R_H) {}

SoundShellSpectrum::SoundShellSpectrum(
    double v_w, 
    double alpha_N, 
    double R_H,
    NucleationHistory nhist,
    ExpansionHistory exhist,
    double step_size, 
    profile_calculator::ProfileSolveMode solve_mode
) : 
    fluid_profile(v_w, alpha_N, step_size, solve_mode), 
    nuc_hist(nhist),
    exp_hist(exhist), 
    H_R(R_H) 
{}

void SoundShellSpectrum::set_chi(double chi_min, double chi_max, size_t N) {
    set_chi(utils::log_space(chi_min, chi_max, N));
}

void SoundShellSpectrum::set_chi(const std::vector<double>& chi_vec) {
    chi = chi_vec;
    A2.clear(); 
    Ekin_over_R.clear(); 
    k_EPi.clear();
    Uf2 = -1.0;
    script_C = -1.0;
    if ( calculate_duration ) { H_delta_tau_fin = -1.0; }
    OmegaGW.clear();
    Delta_tilde_0 = -1.0;
}

void SoundShellSpectrum::set_qR(double qR_min, double qR_max, size_t N) {
    set_qR(utils::log_space(qR_min, qR_max, N));
}

void SoundShellSpectrum::set_qR(const std::vector<double>& qR_vec) {
    qR = qR_vec;
    Ekin_over_R.clear(); 
    k_EPi.clear();
    Uf2 = -1.0;
    script_C = -1.0;
    if ( calculate_duration ) { H_delta_tau_fin = -1.0; }
    OmegaGW.clear();
    Delta_tilde_0 = -1.0;
}

void SoundShellSpectrum::set_qR_from_kR(size_t N) {
    set_qR(utils::log_space(
        (kR.front() * kR.front() / kR[1]) 
            * 0.5 * (1.0 - constants::cs) / constants::cs, 
        (kR.back() * kR.back() / kR[kR.size()-2]) 
            * 0.5 * (1.0 + constants::cs) / constants::cs, 
        N
    ));
}

void SoundShellSpectrum::set_kR(
    double kR_min, double kR_max, size_t N, size_t N_qR
) {
    set_kR(utils::log_space(kR_min, kR_max, N), N_qR);
}

void SoundShellSpectrum::set_kR(
    const std::vector<double>& kR_vec, size_t N_qR
) {
    kR = kR_vec;
    k_EPi.clear();
    OmegaGW.clear();
    if ( N_qR != static_cast<size_t>(-1) ) { set_qR_from_kR(N_qR); }
}

std::pair<std::vector<double>,std::vector<double>> 
SoundShellSpectrum::calculate_A_squared(
    double chi_min, double chi_max, size_t N
) {
    set_chi(chi_min, chi_max, N);
    return calculate_A_squared();
}

std::pair<std::vector<double>,std::vector<double>> 
SoundShellSpectrum::calculate_A_squared(const std::vector<double>& chi_vec) {
    set_chi(chi_vec);
    return calculate_A_squared();
}

std::pair<std::vector<double>,std::vector<double>> 
SoundShellSpectrum::calculate_A_squared() {
    if ( chi.empty() ) { throw std::invalid_argument(
        "The vector chi at which the shape function is calculated is empty."
        " You need to call set_chi first."
    );}
    if ( verbose ) {
        std::cout << "Calculating sine transforms at " << chi.size() 
                  << " points between " << chi.front() << " and " << chi.back() 
                  << "." << std::endl;
    }
    std::vector<double> df, l;
    std::tie(A2, df, l) = sound_shell_model::calculate_sine_transforms(
        chi, 
        fluid_profile.get_xi_coordinate(), 
        fluid_profile.get_fluid_velocity(),
        fluid_profile.get_enthalpy_density(), 
        fluid_profile.get_energy_density(), 
        fluid_profile.get_wall_index(), 
        fluid_profile.get_shock_index()
    );
    return {df, l};
}

void SoundShellSpectrum::calculate_kinetic_spectrum(
    double qR_min, double qR_max, size_t N
) {
    set_qR(qR_min, qR_max, N);
    calculate_kinetic_spectrum();
}

void SoundShellSpectrum::calculate_kinetic_spectrum(
    const std::vector<double>& qR_vec
) {
    set_qR(qR_vec);
    calculate_kinetic_spectrum();
}

void SoundShellSpectrum::calculate_kinetic_spectrum() {
    if ( A2.empty() ) { calculate_A_squared(); }
    if ( qR.empty() ) { throw std::invalid_argument(
        "The dimensionaless wave-number vector qR at which the kinetic spectrum"
        " is calculated is empty. Call set_qR to set the momenta."
    );}

    if ( verbose ) { 
        std::cout << "Calculating kinetic spectrum for "
                  << ( 
                     nuc_hist  == NucleationHistory::EXPONENTIAL ? 
                     "exponential" : "simultaneous" 
                  ) << " nucleation at " << qR.size() << " points between " 
                  << qR.front() << " and " << qR.back() << "."  << std::endl;
    }
    const double vw = fluid_profile.get_wall_velocity();
    Ekin_over_R = sound_shell_model::scaled_kinetic_spectrum(
        qR, chi, A2, vw, nuc_hist
    );
}

void SoundShellSpectrum::calculate_mean_squared_fluid_velocity() {
    Uf2 = sound_shell_model::mean_squared_fluid_velocity(
        qR, get_kinetic_spectrum()
    );
}

void SoundShellSpectrum::calculate_anisotropic_stress_UETC_normalization() {
    script_C = sound_shell_model::anisotropic_stress_power_spectrum_amplitude(
        qR, get_normalized_kinetic_spectrum()
    );
}

void SoundShellSpectrum::calculate_gw_spectrum(
    double kR_min, double kR_max, size_t N, size_t N_qR
) {
    set_kR(kR_min, kR_max, N, N_qR);
    calculate_gw_spectrum();
}

void SoundShellSpectrum::calculate_gw_spectrum(
    const std::vector<double>& kR_vec, size_t N_qR
) {
    set_kR(kR_vec, N_qR);
    calculate_gw_spectrum();
}

void SoundShellSpectrum::calculate_gw_spectrum() {
    if ( Ekin_over_R.empty() ) { calculate_kinetic_spectrum(); }
    if ( kR.empty() ) { throw std::invalid_argument(
        "The dimensionaless wave-number vector kR at which the GW spectrum is"
        " calculated is empty. Call set_kR to set the momenta."
    );}

    std::vector<double> log_qR(qR.size());
    #pragma omp simd
    for (size_t i=0; i<log_qR.size(); ++i) { log_qR[i] = std::log(qR[i]); }
    utils::UniformLinearInterpolator Ekin_over_R_interpolator(
        log_qR, Ekin_over_R
    );

    if ( H_delta_tau_fin <= 0.0 ) { calculate_source_duration(); }
    double Gamma = fluid_profile.get_adiabatic_index();

    if ( verbose ) {
        std::cout << "Calculating GW spectral density at " << kR.size() 
                  << " points between " << kR.front() << " and " << kR.back() 
                  << " for " << ( 
                     nuc_hist == NucleationHistory::EXPONENTIAL ? 
                     "exponential" : "simultaneous" 
                  ) << " nucleation in a " << ( 
                     exp_hist == ExpansionHistory::RADIATION_DOMINATION ? 
                     "radiation-dominated" : "static"
                  ) << " universe, assuming a mean bubble separation of " 
                  << H_R << ", a source duration of " << H_delta_tau_fin 
                  << " in units of the conformal Hubble radius/time, and an"
                  << " adiabatic index of " << Gamma << "." << std::endl;
    }
    OmegaGW = sound_shell_model::gravitational_wave_power_spectrum(
        kR, 
        Ekin_over_R_interpolator, 
        H_delta_tau_fin, 
        H_R, 
        exp_hist, 
        integration, 
        Gamma, 
        N_p, 
        N_z
    );
}

void SoundShellSpectrum::calculate_anisotropic_stress_autocorrelator() {
    if ( Ekin_over_R.empty() ) { calculate_kinetic_spectrum(); }
    if ( kR.empty() ) { throw std::invalid_argument(
        "The dimensionaless wave-number vector kR at which the GW spectrum is"
        " calculated is empty. Call set_kR to set the momenta."
    );}

    std::vector<double> log_qR(qR.size());
    #pragma omp simd
    for (size_t i=0; i<log_qR.size(); ++i) { log_qR[i] = std::log(qR[i]); }
    utils::UniformLinearInterpolator Ekin_over_R_interpolator(
        log_qR, Ekin_over_R
    );
    double w_bar = fluid_profile.get_mean_enthalpy_density();

    if ( verbose ) {
        std::cout << "Calculating anisotropic stress power spectrum at " 
                  << kR.size() << " points between " << kR.front() << " and " 
                  << kR.back() << " for " << ( 
                      nuc_hist == NucleationHistory::EXPONENTIAL ? 
                      "exponential" : "simultaneous"
                  ) << " nucleation with w_bar = " << w_bar << "." << std::endl;
    }

    k_EPi = sound_shell_model::anisotropic_stress_power_spectrum(
        kR, Ekin_over_R_interpolator, w_bar, N_p, N_z
    );
}

void SoundShellSpectrum::calculate_weighted_average_duration_dependence() {
     Delta_tilde_0 = sound_shell_model::weighted_average_Delta_0(
        qR, get_kinetic_spectrum(), get_source_duration(), H_R, exp_hist
    );
}

const FluidProfile& SoundShellSpectrum::get_fluid_profile() const { 
    return fluid_profile; 
}
const std::vector<double>& SoundShellSpectrum::get_chi() const { return chi; }
const std::vector<double>& SoundShellSpectrum::get_qR() const { return qR; }
const std::vector<double>& SoundShellSpectrum::get_kR() const { return kR; }

const std::vector<double>& SoundShellSpectrum::get_shape_function() { 
    if ( A2.empty() ) { calculate_A_squared(); }
    return A2; 
}

const std::vector<double>& SoundShellSpectrum::get_kinetic_spectrum() { 
    if ( Ekin_over_R.empty() ) { calculate_kinetic_spectrum(); }
    return Ekin_over_R; 
}

std::vector<double> SoundShellSpectrum::get_normalized_kinetic_spectrum() {
    return sound_shell_model::kinetic_spectral_shape(get_kinetic_spectrum());
}

std::vector<double> 
SoundShellSpectrum::get_plane_wave_velocity_spectral_density() {
    return utils::rescale_spectrum(
        get_kinetic_spectrum(), qR, -2, 2.0*constants::Pi*constants::Pi
    ); 
}

std::vector<double> SoundShellSpectrum::get_plane_wave_velocity_power_spectrum()
{
    return utils::rescale_spectrum(get_kinetic_spectrum(), qR, 1, 1.0); 
}

std::vector<double> SoundShellSpectrum::get_velocity_spectral_density() { 
    return utils::rescale_spectrum(
        get_kinetic_spectrum(), qR, -2, 4.0*constants::Pi*constants::Pi
    ); 
}

std::vector<double> SoundShellSpectrum::get_velocity_power_spectrum() {
    return utils::rescale_spectrum(get_kinetic_spectrum(), qR, 1, 2.0); 
}

const std::vector<double>& SoundShellSpectrum::get_gravitational_wave_spectrum()
{ 
    if ( OmegaGW.empty() ) { calculate_gw_spectrum(); }
    return OmegaGW;
}

const std::vector<double>& 
SoundShellSpectrum::get_gravitational_wave_power_spectrum() { 
    return get_gravitational_wave_spectrum();
}

std::vector<double> 
SoundShellSpectrum::get_gravitational_wave_spectral_density() {
    return utils::rescale_spectrum(
        get_gravitational_wave_power_spectrum(), 
        kR, 
        -3, 
        2.0*constants::Pi*constants::Pi
    ); 
}

const std::vector<double>& 
SoundShellSpectrum::get_anisotropic_stress_power_spectrum() {
    if ( k_EPi.empty() ) { calculate_anisotropic_stress_autocorrelator(); }
    return k_EPi;
}

std::vector<double> 
SoundShellSpectrum::get_normalized_anisotropic_stress_power_spectrum() {
    return sound_shell_model::anisotropic_stress_power_spectral_shape(
        get_anisotropic_stress_power_spectrum(),
        get_anisotropic_stress_UETC_normalization(), 
        get_kinetic_spectrum_amplitude(),
        fluid_profile.get_mean_enthalpy_density()
    );
}

std::vector<double> 
SoundShellSpectrum::get_anisotropic_stress_spectral_density() {
    return utils::rescale_spectrum(
        get_anisotropic_stress_power_spectrum(),
         kR, 
         -3, 
         2.0*constants::Pi*constants::Pi
    ); 
}

std::vector<double> 
SoundShellSpectrum::get_normalized_gravitational_wave_spectrum() {
    return sound_shell_model::normalized_GW_spectrum(
        get_gravitational_wave_spectrum(),
        get_anisotropic_stress_UETC_normalization(),
        get_kinetic_spectrum_amplitude(),
        get_weighted_average_duration_dependence(),
        fluid_profile.get_adiabatic_index()
    );
}

std::vector<double> SoundShellSpectrum::get_spectral_modification_function() {
    return sound_shell_model::spectral_modification_function(
        get_normalized_gravitational_wave_spectrum(),
        get_normalized_anisotropic_stress_power_spectrum()
    );
}

double SoundShellSpectrum::get_mean_squared_fluid_velocity() {
    if ( Uf2 <= 0.0 ) { calculate_mean_squared_fluid_velocity(); }
    return Uf2;
}

double SoundShellSpectrum::get_kinetic_energy_density_parameter() {
    return sound_shell_model::kinetic_energy_density_parameter(
        get_mean_squared_fluid_velocity(),
        fluid_profile.get_adiabatic_index()
    );
}

double SoundShellSpectrum::get_kinetic_spectrum_amplitude() {
    if ( Ekin_over_R.empty() ) { calculate_kinetic_spectrum(); }
    return *std::max_element(Ekin_over_R.begin(), Ekin_over_R.end());;
}

double SoundShellSpectrum::get_kinetic_spectrum_peak_position() {
    return utils::find_peak(qR, get_kinetic_spectrum()).first;
}

double SoundShellSpectrum::get_normalized_kinetic_density_parameter() {
    return (
        0.5 * get_mean_squared_fluid_velocity() / 
        get_kinetic_spectrum_amplitude()
    );
}

double SoundShellSpectrum::get_anisotropic_stress_UETC_normalization() {
    if ( script_C < 0.0 ) { calculate_anisotropic_stress_UETC_normalization(); }
    return script_C;
}

double SoundShellSpectrum::get_normalized_anisotropic_stress_UETC_amplitude() {
    return utils::find_peak(
        kR, 
        get_normalized_anisotropic_stress_power_spectrum()
    ).second;
}

double 
SoundShellSpectrum::get_normalized_anisotropic_stress_UETC_peak_position() {
    return utils::find_peak(kR, get_anisotropic_stress_power_spectrum()).first;
}

double SoundShellSpectrum::get_weighted_average_duration_dependence() {
    if ( Delta_tilde_0 <= 0.0) { 
        calculate_weighted_average_duration_dependence(); 
    }
    return Delta_tilde_0;
}

void SoundShellSpectrum::set_nucleation_history(NucleationHistory nhist) { 
    nuc_hist = nhist; 
    Ekin_over_R.clear();
    k_EPi.clear();
    Uf2 = -1.0;
    script_C = -1.0; 
    if ( calculate_duration ) { H_delta_tau_fin = -1.0; }
    OmegaGW.clear();
    Delta_tilde_0 = -1.0;
}

NucleationHistory SoundShellSpectrum::get_nucleation_history() const { 
    return nuc_hist; 
}

void SoundShellSpectrum::set_verbose(bool b) { verbose = b; };
bool SoundShellSpectrum::is_verbose() { return verbose; };

void SoundShellSpectrum::recalculate_profile(
    double v_w, 
    double alpha_N, 
    double step_size, 
    profile_calculator::ProfileSolveMode solve_mode
) {
    fluid_profile = FluidProfile(v_w, alpha_N, step_size, solve_mode);
    A2.clear();
    Ekin_over_R.clear();
    k_EPi.clear();
    Uf2 = -1.0;
    script_C = -1.0;
    if ( calculate_duration ) { H_delta_tau_fin = -1.0; }
    OmegaGW.clear();
    Delta_tilde_0 = -1.0;
}

double SoundShellSpectrum::get_wall_velocity() const { 
    return fluid_profile.get_wall_velocity(); 
}
double SoundShellSpectrum::get_transition_strength() const { 
    return fluid_profile.get_transition_strength(); 
}

void SoundShellSpectrum::set_expansion_history(ExpansionHistory ehist) {
    if ( ehist == ExpansionHistory::NONE ) { throw std::invalid_argument(
        "'NONE' is not a valid expansion history. Please choose"
        " RADIATION_DOMINATION or STATIC_UNIVERSE."
    ); }
    exp_hist = ehist;
    OmegaGW.clear();
    Delta_tilde_0 = -1.0;
}
ExpansionHistory SoundShellSpectrum::get_expansion_history() const { 
    return exp_hist; 
}

void SoundShellSpectrum::set_mean_bubble_separation(double R_H) { 
    H_R = R_H; 
    if ( calculate_duration ) { H_delta_tau_fin = -1.0; }
    OmegaGW.clear();
    Delta_tilde_0 = -1.0;
}

double SoundShellSpectrum::get_mean_bubble_separation() const { return H_R; }

void SoundShellSpectrum::set_source_duration(double duration) {
    calculate_duration = false;
    H_delta_tau_fin = duration;
    OmegaGW.clear();
    Delta_tilde_0 = -1.0;
}

void SoundShellSpectrum::calculate_source_duration() {
    calculate_duration = true;
    H_delta_tau_fin = H_R/std::sqrt(
        get_kinetic_energy_density_parameter()
    );
    OmegaGW.clear();
    Delta_tilde_0 = -1.0;
}

double SoundShellSpectrum::get_source_duration() { 
    if ( H_delta_tau_fin <= 0.0 ) { calculate_source_duration(); }
    return H_delta_tau_fin;
}

void SoundShellSpectrum::set_correlator_integration_settings(
    size_t n_p, size_t n_z, sound_shell_model::MomentumIntegration integ
) {
    OmegaGW.clear(); 
    if ( n_p != N_p || n_z != N_z ) { k_EPi.clear(); }
    N_p = n_p; N_z = n_z; integration = integ;
}

} // namespace inverse_pt