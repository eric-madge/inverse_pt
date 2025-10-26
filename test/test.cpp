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

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>

#include "inverse_pt/fluid_profile.hpp"
#include "inverse_pt/sound_shell_spectrum.hpp"
#include "inverse_pt/profile_calculator.hpp"
#include "inverse_pt/sound_shell_model.hpp"
#include "inverse_pt/settings.hpp"
#include "inverse_pt/constants.hpp"
#include "inverse_pt/utils.hpp"

using namespace inverse_pt::utils;
using namespace inverse_pt::profile_calculator;
using namespace inverse_pt::sound_shell_model;
using namespace inverse_pt::constants;

int main(int argc, char* argv[]) {

    //inverse_pt::settings::debug = true;
    //inverse_pt::settings::check_mode = inverse_pt::settings::CheckMode::ERROR;

    if (argc < 3) {
        std::cerr << "Usage: ./test <vw> <alpha_N> [HR] [Hdtau] [precision]\n";
        std::cerr << "Example: ./test 0.4 0.1 0.01 1.0 0.001\n";
        return 1;
    }
    double vw = std::stod(argv[1]);
    double alpha_N = std::stod(argv[2]);
    double H_R = (argc >= 4) ? std::stod(argv[3]) : 1.0;
    double H_delta_tau = (argc >= 5) ? std::stod(argv[4]) : 1.0;
    double precision = (argc >= 6) ? std::stod(argv[5]) : 1e-3;

    std::chrono::time_point<std::chrono::high_resolution_clock> main_start 
        = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> sub_start;
    std::chrono::duration<double> duration;

    std::vector<double> xi, v;

    std::cout << "Calculating fluid profile for xi_w = " << vw 
              << " and alpha_N = " << alpha_N << "." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();

    HydroSolutionType sol_type = classify_wave_fluid_frame(
        vw, alpha_N, precision
    );
    double alpha_plus = find_alpha_plus(sol_type, alpha_N, vw, precision);
    auto [vplus, vminus] = compute_fluid_velocities(vw, alpha_plus, sol_type);

    std::cout << "Found " << sol_type << " solution with alpha_+ = " 
              << alpha_plus << ".\nFluid velocities: v+ = " << vplus 
              << ", v_- = " << vminus << ".\nCalculating velocity profile."
              << std::endl;

    size_t idx_wall = compute_velocity_profile(
        sol_type, xi,v, vw, alpha_plus, precision
    );
    size_t idx_shock = NO_INDEX;
    if ( 
        sol_type != HydroSolutionType::Detonation 
        && sol_type != HydroSolutionType::InverseDetonation 
    ) {
        idx_shock = get_shock_index(sol_type, xi.size());
    }
    std::cout << "Velocity profile calculated with " << xi.size() << " points."
              << "\n  Velocity jumps from " << v[idx_wall-1] << " to " 
              << v[idx_wall] << " at the bubble wall, xi = " << xi[idx_wall] 
              << " (index: " << idx_wall << ")." ;
    if ( idx_shock != NO_INDEX ) { 
        std::cout << "\n  Velocity jumps from " << v[idx_shock-1] << " to " 
                  << v[idx_shock] << " at the shock position, xi = " 
                  << xi[idx_shock] << " (index: " << idx_shock << ").";
    }
    std::cout << std::endl;
    std::vector<double> w = compute_enthalpy_density(
        xi, v, idx_wall, idx_shock
    );
    double wN = w.back();
    std::cout << "Enthalpy density calculated. wN: " << wN << std::endl;

    double epsilon = compute_epsilon(wN, alpha_N);
    std::cout << "Vacuum energy: epsilon = " << epsilon << std::endl;
    std::vector<double> e = compute_energy_density_profile(
        w, idx_wall, epsilon
    );
    std::vector<double> p = compute_pressure_density_profile(
        w, idx_wall, epsilon
    );

    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "All profiles calculated after " << duration.count() 
              << " seconds. Writing to file." << std::endl;
    std::ofstream outfile("profiles.csv");
    outfile << "# xi, v, w, e, p\n";
    for ( size_t i=0; i<xi.size(); ++i ) {
        outfile << xi[i] << ", " << v[i] << ", " << w[i] << ", " << e[i] << ", "
                << p[i] << "\n";
    }
    outfile.close();

    inverse_pt::FluidProfile::set_verbose(true);
    inverse_pt::FluidProfile fluid_profile(xi, v, w, e, p, idx_wall, idx_shock);

    std::cout << "Calculating GW spectrum with H*R = " << H_R 
              << " and H*deltaTau = " << H_delta_tau << ":" << std::endl;

    inverse_pt::SoundShellSpectrum::set_verbose(true);
    inverse_pt::SoundShellSpectrum spectra(fluid_profile, H_R);

    spectra.set_chi(1e-2, 1e4, 2401); 
    spectra.set_kR(1e-3, 1e4, 200);
    spectra.set_qR_from_kR(10001);

    std::vector<double> chi = spectra.get_chi();
    std::vector<double> qR = spectra.get_qR();
    std::vector<double> kR = spectra.get_kR();

    if ( H_delta_tau > 0. ) { spectra.set_source_duration(H_delta_tau); }
    size_t N_p = 200, N_z = 64;
    spectra.set_correlator_integration_settings(N_p, N_z);

    sub_start = std::chrono::high_resolution_clock::now();
    auto [df, l] = spectra.calculate_A_squared();
    std::vector<double> A2 = spectra.get_shape_function();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Shape function calculated after " << duration.count() 
              << " seconds." << std::endl;

    spectra.set_nucleation_history(NucleationHistory::EXPONENTIAL);
    std::cout << "Calculating spectra for exponential nucleation:" << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> E_kin_over_R_exp = spectra.get_kinetic_spectrum();
    double K_peak_kin_exp = spectra.get_kinetic_spectrum_peak_position();
    double E_kin_star_over_R_exp = spectra.get_kinetic_spectrum_amplitude();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Kinetic spectrum calculated after " << duration.count() 
              << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    double Uf2_exp = spectra.get_mean_squared_fluid_velocity();
    double Omega_K_exp = spectra.get_kinetic_energy_density_parameter();
    double curly_K_exp = spectra.get_normalized_kinetic_density_parameter();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Kinetic density parameter calculated after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> k_E_Pi_exp = 
        spectra.get_anisotropic_stress_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Anisotropic stress power spectrum calculated after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    double curly_C_exp = spectra.get_anisotropic_stress_UETC_normalization();
    std::vector<double> K3_zeta_Pi_exp = 
        spectra.get_normalized_anisotropic_stress_power_spectrum();
    auto [K_peak_Pi_exp, K3_zeta_Pi_peak_exp] = find_peak(
        spectra.get_kR(), K3_zeta_Pi_exp
    );
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Normalized anisotropic stress calculated after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_expansion_history(ExpansionHistory::STATIC_UNIVERSE);
    double tilde_Delta_0_flat_exp = 
        spectra.get_weighted_average_duration_dependence();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated amplitude modification in a flat universe after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_correlator_integration_settings(
        N_p, N_z, MomentumIntegration::INFINITE_DURATION
    );
    std::vector<double> Omega_GW_approx_exp = 
        spectra.get_gravitational_wave_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated approximate GW spectrum after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_correlator_integration_settings(
        N_p, N_z, MomentumIntegration::FULL
    );
    std::vector<double> Omega_GW_flat_exp = 
        spectra.get_gravitational_wave_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated GW spectrum in a flat universe after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> K3_zeta_GW_flat_exp = 
        spectra.get_normalized_gravitational_wave_spectrum();
    std::vector<double> tilde_Delta_flat_exp =
        spectra.get_spectral_modification_function();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Normalized GW spectrum in a flat universe after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_expansion_history(ExpansionHistory::RADIATION_DOMINATION);
    double tilde_Delta_0_rad_exp = 
        spectra.get_weighted_average_duration_dependence();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated amplitude modification in an expanding universe"
              << " after " << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> Omega_GW_rad_exp = 
        spectra.get_gravitational_wave_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated GW spectrum in an expanding universe after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> K3_zeta_GW_rad_exp = 
        spectra.get_normalized_gravitational_wave_spectrum();
    std::vector<double> tilde_Delta_rad_exp = 
        spectra.get_spectral_modification_function();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Normalized GW spectrum in an expanding universe after " 
              << duration.count() << " seconds." << std::endl;

    spectra.set_nucleation_history(NucleationHistory::SIMULTANEOUS);
    std::cout << "Calculating spectra for simultaneous nucleation:" 
              << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> E_kin_over_R_sim = spectra.get_kinetic_spectrum();
    double K_peak_kin_sim = spectra.get_kinetic_spectrum_peak_position();
    double E_kin_star_over_R_sim = spectra.get_kinetic_spectrum_amplitude();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Kinetic spectrum calculated after " << duration.count() 
              << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    double Uf2_sim = spectra.get_mean_squared_fluid_velocity();
    double Omega_K_sim = spectra.get_kinetic_energy_density_parameter();
    double curly_K_sim = spectra.get_normalized_kinetic_density_parameter();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Kinetic density parameter calculated after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> k_E_Pi_sim = 
        spectra.get_anisotropic_stress_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Anisotropic stress power spectrum calculated after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    double curly_C_sim = spectra.get_anisotropic_stress_UETC_normalization();
    std::vector<double> K3_zeta_Pi_sim =
        spectra.get_normalized_anisotropic_stress_power_spectrum();
    auto [K_peak_Pi_sim, K3_zeta_Pi_peak_sim] = find_peak(
        spectra.get_kR(), K3_zeta_Pi_sim
    );
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Normalized anisotropic stress calculated after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_expansion_history(ExpansionHistory::STATIC_UNIVERSE);
    double tilde_Delta_0_flat_sim = 
        spectra.get_weighted_average_duration_dependence();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated amplitude modification in a flat universe after "
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_correlator_integration_settings(
        N_p, N_z, MomentumIntegration::INFINITE_DURATION
    );
    std::vector<double> Omega_GW_approx_sim = 
        spectra.get_gravitational_wave_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated approximate GW spectrum after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_correlator_integration_settings(
        N_p, N_z, MomentumIntegration::FULL
    );
    std::vector<double> Omega_GW_flat_sim = 
        spectra.get_gravitational_wave_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated GW spectrum in a flat universe after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> K3_zeta_GW_flat_sim = 
        spectra.get_normalized_gravitational_wave_spectrum();
    std::vector<double> tilde_Delta_flat_sim = 
        spectra.get_spectral_modification_function();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Normalized GW spectrum in a flat universe after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    spectra.set_expansion_history(ExpansionHistory::RADIATION_DOMINATION);
    double tilde_Delta_0_rad_sim = 
        spectra.get_weighted_average_duration_dependence();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated amplitude modification in an expanding universe"
              << " after " << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> Omega_GW_rad_sim = 
        spectra.get_gravitational_wave_power_spectrum();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Calculated GW spectrum in an expanding universe after " 
              << duration.count() << " seconds." << std::endl;
    sub_start = std::chrono::high_resolution_clock::now();
    std::vector<double> K3_zeta_GW_rad_sim = 
        spectra.get_normalized_gravitational_wave_spectrum();
    std::vector<double> tilde_Delta_rad_sim = 
        spectra.get_spectral_modification_function();
    duration = std::chrono::high_resolution_clock::now() - sub_start;
    std::cout << "  Normalized GW spectrum in an expanding universe after "
              << duration.count() << " seconds." << std::endl;

    duration = std::chrono::high_resolution_clock::now() - main_start;
    std::cout << "\nTotal run time: " << duration.count() << " seconds.\n"
              << "Results:\n"
              << "  exponential nuclation:\n"
              << "    10^4 E_kin_*/R_* = " << 1.0e4*E_kin_star_over_R_exp <<"\n"
              << "    K_kin_peak       = " << K_peak_kin_exp << "\n"
              << "    curly_K          = " << curly_K_exp << "\n"
              << "    200 Uf^2         = " << 2.0e2*Uf2_exp << "\n"
              << "    10^2 Omega_K     = " << 1.0e2*Omega_K_exp << "\n"
              << "    curly_C          = " << curly_C_exp << "\n"
              << "    K_GW_peak        = " << K_peak_Pi_exp << "\n"
              << "    K^3 zeta_Pi_peak = " << K3_zeta_Pi_peak_exp << "\n"
              << "    Delta~0 (rad.)   = " << tilde_Delta_0_rad_exp << "\n"
              << "    Delta~0 (stat.)  = " << tilde_Delta_0_flat_exp << "\n"
              << "  simultaneous nuclation:\n"
              << "    10^4 E_kin_*/R_* = " << 1.0e4*E_kin_star_over_R_sim <<"\n"
              << "    K_kin_peak       = " << K_peak_kin_sim << "\n"
              << "    curly_K          = " << curly_K_sim << "\n"
              << "    200 Uf^2         = " << 2.0e2*Uf2_sim << "\n"
              << "    10^2 Omega_K     = " << 1.0e2*Omega_K_sim << "\n"
              << "    curly_C          = " << curly_C_sim << "\n"
              << "    K_GW_peak        = " << K_peak_Pi_sim << "\n"
              << "    K^3 zeta_Pi_peak = " << K3_zeta_Pi_peak_sim << "\n"
              << "    Delta~0 (rad.)   = " << tilde_Delta_0_rad_sim << "\n"
              << "    Delta~0 (stat.)  = " << tilde_Delta_0_flat_sim << "\n"
              << std::endl;

    std::cout << "Done. Writing spectra to file" << std::endl;
    outfile.open("spectra.csv");
    outfile << "# chi, df_dchi, l, A2, qR, Ekin_exp, Ekin_sim, kR, kEPi_exp, "
            << "kEPi_sim, OmegaGW_approx_exp, OmegaGW_approx_sim, "
            << "OmegaGW_flat_exp, OmegaGW_flat_sim, OmegaGW_exp, OmegaGW_sim\n";
    for ( size_t i=0; i<kR.size() || i<qR.size() || i<chi.size(); ++i ) {
        if ( i < chi.size() ) { 
            outfile << chi[i] << ", " << df[i] << ", " << l[i] << ", "
                    << A2[i] << ", "; 
        }
        else { outfile << ",,,, "; }
        if ( i < qR.size() ) { 
            outfile << qR[i] << ", " << E_kin_over_R_exp[i] << ", " 
                    << E_kin_over_R_sim[i] << ", "; 
        }
        else { outfile << ",,, "; }
        if ( i < kR.size() ) { 
            outfile << kR[i] << ", " << k_E_Pi_exp[i]
                             << ", " << k_E_Pi_sim[i] 
                             << ", " << Omega_GW_approx_exp[i] 
                             << ", " << Omega_GW_approx_sim[i] 
                             << ", " << Omega_GW_flat_exp[i] 
                             << ", " << Omega_GW_flat_sim[i] 
                             << ", " << Omega_GW_rad_exp[i] 
                             << ", " << Omega_GW_rad_sim[i] 
                             << "\n"; 
        } else { outfile << ",,,,,,,,\n"; }
    }
    outfile.close();
    
    std::cout << "Done." << std::endl;

    return 0;
}