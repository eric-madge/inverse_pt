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

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "inverse_pt_calculator_bindings.hpp"
#include "python_helpers.hpp"
#include "inverse_pt/sound_shell_model.hpp"

namespace py = pybind11;
using namespace inverse_pt;
using namespace inverse_pt::sound_shell_model;
using MomentumIntegration 
    = AnisotropicStressCorrelatorSettings::MomentumIntegration;

void bind_sound_shell_model(pybind11::module_& m) {

    // Nucleation history types
    py::enum_<NucleationHistory>(
        m, "NucleationHistory", py::arithmetic(), "Nucleation history options"
    )
        .value(
            "EXPONENTIAL",  
            NucleationHistory::EXPONENTIAL,  
            "exponential nucleation"
        )
        .value(
            "SIMULTANEOUS", 
            NucleationHistory::SIMULTANEOUS, 
            "simultaneous nucleation"
        )
    .export_values();

    // Expansion history types
    py::enum_<ExpansionHistory>(
        m, "ExpansionHistory", py::arithmetic(), "Expansion history options"
    )
        .value(
            "RADIATION_DOMINATION", 
            ExpansionHistory::RADIATION_DOMINATION, 
            "assume a radiation dominated universe"
        )
        .value(
            "STATIC_UNIVERSE", 
            ExpansionHistory::STATIC_UNIVERSE,
            "assume a static (non-expanding) universe"
        )
        .value(
            "NONE", 
            ExpansionHistory::NONE, 
            "no expansion history / time integration (i.e. use Delta=1 to"
            " calculate the autocorrelation zeta_Pi"
        )
    .export_values();

    // Approximation used in the momentum integration
    py::enum_<MomentumIntegration>(
        m, 
        "MomentumIntegration", 
        py::arithmetic(), 
        "Approximation used in the momentum integration of the anisotropic"
        " stress corrleator"
    )
        .value(
            "FULL", 
            MomentumIntegration::FULL, 
            "Full momentum integration, using Eqs. (56) and (59) of [Roper+"
            " (2024)] for Delta and Delta_flat, respectively."
        )
        .value(
            "LOW_FREQUENCY", 
            MomentumIntegration::LOW_FREQUENCY, 
            "Use the low-frequency limit Delta_0 or Delta_0_flat, i.e. Eqs."
            " (61,62) of [Roper+ (2024)] for Delta and Delta_flat."
        )
        .value(
            "INFINITE_DURATION",  
            MomentumIntegration::INFINITE_DURATION, 
            "Use the inifinite-duration (high-frequency) limit of Delta,"
            " k * delta_tau_fin >> 1, i.e. Eqs. (86) and (83) of [Roper+"
            " (2024)], respectively. "
        )
    .export_values();

    py::class_<AnisotropicStressCorrelatorSettings> (
        m, 
        "AnisotropicStressCorrelatorSettings",
        "Settings for the momentum and time integration of the anisotropic"
        " stress and gravitational wave spectrum."
    )
        .def(
            py::init<>(),
            "Default constructor (for calculating the anisotropic stress"
            " autocorrelation power spectrum k E_Pi(k)."
        )
        .def(
            py::init<size_t, size_t, double>(),
            py::arg("N_p"), py::arg("N_z"), py::arg("c_w")=2.0,
            "Constuctor for settings for calculating the anisotropic stress"
            " autocorrelation power spectrum k E_Pi(k).\n\n"
            "Note: This sets universe to ExpansionHistory.NONE and freq_approx"
            " to MomentumIntegration.FULL.\n"
            "\nParameters\n----------\n"
            "N_p : integer\n    Number of points in the integration over"
            " P = p R*\n"
            "N_z : integer\n    Number of points in the integration over z\n"
            "c_w : float, optional (default: 2.0)\n"
            "    Numeric prefactor of the integral in anisotropic stress"
            " correlator calculations."
        )
        .def(
            py::init<
                double, 
                double, 
                ExpansionHistory, 
                MomentumIntegration, 
                double
            >(),
            py::arg("H_delta_tau"), 
            py::arg("H_R"),
            py::arg("universe")=ExpansionHistory::RADIATION_DOMINATION,
            py::arg("freq_approx")=MomentumIntegration::FULL,
            py::arg("c_w")=16.0/3.0,
            "Constructor for time integration settings for the graviational"
            " wave spectrum Omega_GW.\n\n"
            "Note: This automatically sets c_w to 16.0/3.0.\n"
            "\nParameters\n----------\n"
            "H_delta_tau : float\n    Source duration in units of the Hubble"
            " time, H* delta_tau_fin\n"
            "H_R : float\n    Mean bubble separation in units of the Hubble"
            " size, H* R*\n"
            "universe : ExpansionHistory, optional (default:"
            " ExpansionHistory.RADIATION_DOMINATION)\n"
            "    Assumptions used in the time integration.\n"
            "freq_approx : MomentumIntegration, optional (default:"
            " MomentumIntegration.FULL)\n"
            "    Approximation used in the momentum integration of the time"
            " integrals.\n"
            "c_w : float, optional (default: 16.0/3.0)\n"
            "    Numeric prefactor of the integral in anisotropic stress"
            " correlator calculations."
        )
        .def(
            py::init<
                size_t, 
                size_t, 
                double, 
                double, 
                ExpansionHistory,
                MomentumIntegration,
                double
            >(),
            py::arg("N_p"), 
            py::arg("N_z"),
            py::arg("H_delta_tau"), 
            py::arg("H_R"),
            py::arg("universe")=ExpansionHistory::RADIATION_DOMINATION,
            py::arg("freq_approx")=MomentumIntegration::FULL,
            py::arg("c_w")=16.0/3.0,
            "Full constructor.\n"
            "\nParameters\n----------\n"
            "N_p : integer\n"
            "    Number of points in the integration over P = p R*\n"
            "N_z : integer\n"
            "    Number of points in the integration over z\n"
            "H_delta_tau : float\n    Source duration in units of the Hubble"
            " time, H* delta_tau_fin\n"
            "H_R : float\n"
            "    Mean bubble separation in units of the Hubble size, H* R*\n"
            "universe : ExpansionHistory, optional (default: "
            "ExpansionHistory.RADIATION_DOMINATION)\n"
            "    Assumptions used in the time integration.\n"
            "freq_approx : MomentumIntegration, optional"" (default:"
            " MomentumIntegration.FULL)\n    Approximation used in the momentum"
            " integration of the time integrals.\n"
            "c_w : float, optional (default: 16.0/3.0)\n    Numeric prefactor"
            " of the integral in anisotropic stress correlator calculations."
        )
        .def_readwrite(
            "N_p", 
            &AnisotropicStressCorrelatorSettings::N_p, 
            "Number of points in the integration over p R_* (default: 5000)"
        )
        .def_readwrite(
            "N_z", 
            &AnisotropicStressCorrelatorSettings::N_z, 
            "Number of points in the integration over z (default: 100)"
        )
        .def_readwrite(
            "H_delta_tau", 
            &AnisotropicStressCorrelatorSettings::H_delta_tau, 
            "Source duration in units of the Hubble time, H* delta_tau_fin"
        )
        .def_readwrite(
            "H_R", 
            &AnisotropicStressCorrelatorSettings::H_R, 
            "Mean bubble separation in units of the Hubble size, H* R*"
        )
        .def_readwrite(
            "universe", 
            &AnisotropicStressCorrelatorSettings::universe, 
            "Assumptions used in the time integration "
        )
        .def_readwrite(
            "freq_approx", 
            &AnisotropicStressCorrelatorSettings::freq_approx,
            "Approximation used in the momentum integration of the time"
            " integrals"
        )
        .def_readwrite(
            "c_w", 
            &AnisotropicStressCorrelatorSettings::c_w, 
            "Numeric prefactor of the integral in"
            " anisotropic_stress_correlator_integration"
        )
    ;  

    m.def(
        "nu_exp",
        &nu_exp,
        py::arg("T"),
        "Bubble lifetime distribution function for exponential nucleation.\n\n"
        "This is Eq. (4.27) of [Hindmarsh+ (2019)].\n"
        "\nParameters\n----------\n"
        "T : float\n    Dimensionless bubble lifetime T~ = beta T_i\n"
        "\nReturns\n-------\n"
        "nu_exp : float\n    Bubble lifetime distribution evaluated at T for"
        " exponential nucleation"
    );

    m.def(
        "nu_sim",
        &nu_sim,
        py::arg("T"),
        "Bubble lifetime distribution function for simultaneous nucleation.\n\n"
        "This is Eq. (4.32) of [Hindmarsh+ (2019)].\n"
        "\nParameters\n----------\n"
        "T : float\n    Dimensionless bubble lifetime T~ = beta T_i\n"
        "\nReturns\n-------\n"
        "nu_exp : float\n    Bubble lifetime distribution evaluated at T for"
        " simultaneous nucleation"
    );

    m.def(
        "calculate_sine_transforms",
        [](
            py::array_t<double> chi_arr,
            py::array_t<double> xi,
            py::array_t<double> v,
            py::array_t<double> w,
            py::array_t<double> e,
            size_t idx_wall,
            py::object idx_sh
        ) {
            std::vector<double> A2, df, l;
            std::tie(A2, df, l) = calculate_sine_transforms(
                np2vec(chi_arr), 
                np2vec(xi), 
                np2vec(v), 
                np2vec(w), 
                np2vec(e), 
                idx_wall, 
                obj2idx(idx_sh)
            );
            return py::make_tuple(vec2np(A2), vec2np(df), vec2np(l));
        },
        py::arg("chi_arr"),
        py::arg("xi"), 
        py::arg("v"), 
        py::arg("w"), 
        py::arg("e"), 
        py::arg("idx_wall"), 
        py::arg("idx_sh")=py::none(),
        "Calculate the (derivative of) sine transforms of the velocity and"
        " energy spectrum.\n\n"
        "This function calculates the derivative of the sine transform f(chi)"
        " of the velocity spectrum v, the sine transform l(chi) of the energy"
        " fluctuation variable lambda = (e-e_bar)/w_bar, and the absolute value"
        " squared of their sum, |A(chi)|^2 = 1/4 [(f'(chi))^2 + (cs l(chi))^2],"
        " given in Eqs. (29-31) of [Roper+ (2024)] (or Eqs. (4.5), (4.11) and"
        " (4.8) of [Hindmarsh+ (2019)]).\n"
        "\nParameters\n----------\n"
        "chi_arr : numpy.ndarray, shape (N,)\n"
        "    Array of values of the dimensionless wavenumber chi = k Tn\n"
        "xi : numpy.ndarray, shape (M,)\n"
        "    Values of the self-similar coordinate xi\n"
        "v :  numpy.ndarray, shape (M,)\n    Fluid velocity profile v(xi)\n"
        "w :  numpy.ndarray, shape (M,)\n    Enthalpy density profile w(xi)\n"
        "e :  numpy.ndarray, shape (M,)\n    Energy density profile e(xi)\n"
        "idx_wall : integer\n""    Index of the bubble wall position"
        " (first element outside the wall)\n"
        "idx_sh : integer, optional (default: None)\n"
        "    Index of the shock front position (first element outside the"
        " shock, use None if there is no shock front)\n"
        "\nReturns\n-------\n"
        "A2 : numpy.ndarray, shape (N,)\n    Shape function |A(chi)|\n"
        "df : numpy.ndarray, shape (N,)\n    Derivative f'(chi) of the sine"
        " transform of the velocity spectrum\n"
        "l : numpy.ndarray, shape (N,)\n"
        "    Sine transform l(chi) of the energy fluctuation variable"
    );

    m.def(
        "scaled_kinetic_spectrum",
        [](
            py::array_t<double> kR, 
            py::array_t<double> chi_arr, 
            py::array_t<double> A2_arr, 
            double xi_wall, 
            NucleationHistory nuc_hist
        ) {
            return vec2np(scaled_kinetic_spectrum(
                np2vec(kR), np2vec(chi_arr), np2vec(A2_arr), xi_wall, nuc_hist
            ));
        },
        py::arg("kR"), 
        py::arg("chi_arr"), 
        py::arg("A2_arr"), 
        py::arg("xi_wall"), 
        py::arg("nuc_hist")=NucleationHistory::EXPONENTIAL,
        "Scaled kinetic spectrum Ekin(k)/R*.\n\n"
        "This is the momentum dependent coefficient Ekin(k) of the velocity"
        " field unequal-time correlator,\n"
        "    Ekin(k, tau1, tau2) = Ekin(k) cos(k cs (tau1-tau2))\n"
        "cf. Eqs. (39), (33) and (35) of [Roper+ (2024)]. We use the"
        " dimensionless momentum z = k R* as input variable and divide by $R_*$"
        " to make the kinetic spectrum dimensionless.\n\n"
        "Note that this is related to the velocity spectral density of the"
        " plane wave components of the velocity field, Eq. (4.17) of"
        " [Hindmarsh+ (2019)], S_v(k) = 2 pi^2 Ekin(k)/k^2.\n"
        "\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n    Momentum k times bubble separation"
        " R* at which the kinetic spectrum is evaluated\n"
        "chi_arr : numpy.ndarray, shape (M,)\n"
        "    Values of chi = T q = T~ q/beta at which |A|^2 is provided\n"
        "A2_arr : numpy.ndarray, shape (M,)\n    Precomputed values of the"
        " shape function |A|^2 used in the integration\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "nuc_hist : NucleationHistory, optional (default:"
        " NucleationHistory.EXPONENTIAL)\n    Nucleation history (exponential"
        " of simultaneous) used for the collision time distribution\n"
        "\nReturns\n-------\n"
        "Ekin_over_R : numpy.ndarray, shape (N,)\n"
        "    Scaled kinetic spectrum Ekin(k R_*)/R*"
    );

    m.def(
        "kinetic_spectral_shape",
        [](py::array_t<double> E_kin_over_R) { 
            return vec2np(kinetic_spectral_shape(np2vec(E_kin_over_R))); 
        },
        py::arg("E_kin_over_R"),
        "Spectral shape zeta_kin(k R*) of the kinetic spectrum.\n\n"
        "This is the spectral shape zeta_kin(k R*) of the kinetic spectrum,"
        " Ekin(k) = Ekin* zeta_kin(k R*) (cf. Eq. (41) of [Roper+ (2024)],"
        " where Ekin* is the maximum amplitude of the  kinetic spectrum and R*"
        " is the mean bubble separation.\n"
        "\nParameters\n----------\n"
        "E_kin_over_R : numpy.ndarray, shape (N,)\n"
        "    Scaled kinetic spectrum Ekin(k)/R*\n"
        "\nReturns\n-------\n"
        "zeta_kin : numpy.ndarray, shape (N,)\n"
        "    kinetic spectral shape zeta_kin(k R*) = Ekin(k)/Ekin*"
    );

    m.def(
        "mean_squared_fluid_velocity",
        [](py::array_t<double> kR, py::array_t<double> E_kin_over_R) { 
            return mean_squared_fluid_velocity(
                np2vec(kR), np2vec(E_kin_over_R)
            ); 
        },
        py::arg("kR"), py::arg("E_kin_over_R"),
        "Mean square fluid velocity squared Uf^2.\n\n"
        "This is the root-mean-square fluid velocity squared, cf. Eqs. (4.33)"
        " and (4.34) of [Hindmarsh+ (2019)], and footnote 6 of [Roper+ (2024)]," 
        " calculated from the kinetic spectrum (or velocity-field UETC at equal"
        " times, see also Eq. (42) of [Roper+ (2024)]).\n"
        "\nNote\n----\nThe sound shell model assumes non-relativistic fluid"
        " velocities and a constant entropy density, w(x) = w_bar ~ wN, where"
        " wN is the enthalphy density at the nucleation temperature. Thus, the"
        " the root-mean-square fluid velocity coincides with the"
        " enthalpy-weighted average fluid velocity.\n"
        "\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n    Momentum k times bubble separation"
        " R* at which the kinetic spectrum is evaluated\n"
        "E_kin_over_R : numpy.ndarray, shape (N,)\n    Scaled kinetic spectrum"
        " Ekin(k)/R*\n"
        "\nReturns\n-------\n"
        "Uf2 : float\n    Mean-square fluid velocity Uf^2"
    );

    m.def(
        "mean_squared_fluid_velocity",
        [](
            py::array_t<double> chi_arr, 
            py::array_t<double> A2_arr, 
            double xi_wall
        ) { 
            return mean_squared_fluid_velocity(
                np2vec(chi_arr), np2vec(A2_arr), xi_wall
            ); 
        },
        py::arg("chi_arr"), py::arg("A2_arr"), py::arg("xi_wall"),
        "Mean square fluid velocity squared Uf^2.\n\n"
        "This is the root-mean-square fluid velocity squared, cf. Eqs. (4.33)"
        " and (4.34) of [Hindmarsh+ (2019)].\n"
        "The integral over T~ in Eq. (4.33) has been carried out analytically,"
        " assuming that the nucleation history is either exponential or"
        " simultaneous [Hindmarsh+ (2019)].\n"
        "\nParameters\n----------\n"
        "chi_arr : numpy.ndarray, shape (N,)\n"
        "    Values of chi = T q = T~ q/beta at which |A|^2 is provided\n"
        "A2_arr : numpy.ndarray, shape (N,)\n    Precomputed values of the"
        " shape function |A|^2 used in the integration\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "\nReturns\n-------\n"
        "Uf2 : float\n    Mean-square fluid velocity Uf^2"
    );

    m.def(
        "kinetic_energy_density_parameter",
        py::overload_cast<double,double>(&kinetic_energy_density_parameter),
        py::arg("Uf2"), py::arg("Gamma")=4.0/3.0,
        "Total kinetic energy density parameter OmegaK.\n\n"
        "This is the total kinetic energy density OmegaK = rho_kin/rho_crit as"
        " a fraction of the critical energy density, where the kinetic energy"
        " is given by rho_kin = <w gamma^2 v^2>. It is related to the"
        " root-mean-square fluid velocity by OmegaK = Gammma Uf^2,  where"
        " Gamma = w_bar/e_bar is the mean adiabatic index of the fluid in the"
        " stable phase.\n"
        "\nNote\n----\n"
        "This assumes Gamma=4/3 by default. To obtain OmegaK = 1/2 Uf^2 as"
        " defined in Eq. (42) of [Roper+ (2024)], set Gamma to 0.5.\n"
        "\nParameters\n----------\n"
        "Uf2 : float\n    Mean-square fluid velocity Uf^2"
        "Gamma : float, optional (default: 4.0/3.0)\n    Mean adiabatic index"
        " Gamma = w_bar/e_bar of the fluid in the stable phase\n"
        "\nReturns\n-------\n"
        "OmegaK : float\n    kinetic energy density parameter OmegaK\n"
    );

    m.def(
        "kinetic_energy_density_parameter",
        [](
            py::array_t<double> kR, 
            py::array_t<double> E_kin_over_R, 
            double Gamma
        ) { 
            return kinetic_energy_density_parameter(
                np2vec(kR), np2vec(E_kin_over_R), Gamma
            );
        },
        py::arg("kR"), py::arg("E_kin_over_R"), py::arg("Gamma")=4.0/3.0,
        "Total kinetic energy density parameter OmegaK.\n\n"
        "This is the total kinetic energy density OmegaK = rho_kin/rho_crit as"
        " a fraction of the critical energy density, calculated from the"
        " kinetic spectrum (or velocity-field UETC at equal times, cf. Eq. (42)"
        " of [Roper+ (2024)]).\n"
        "\nNote\n----\n-) Set Gamma to 0.5 to obtain OmegaK as defined in Eq."
        " (42) of [Roper+ (2024)].\n"
        "-) To obtain the integral of the kinetic spectral shape defined in Eq."
        " (43) of [Roper+ (2024)], script_K, provide the kinetic spectral shape"
        " zeta_lin(k R*) = Ekin(k)/Ekin* for E_kin_over_R and set Gamma to 0.5."
        "\n\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n    Momentum k times bubble separation"
        " R* at which the kinetic spectrum is evaluated\n"
        "E_kin_over_R : numpy.ndarray, shape (N,)\n"
        "    Scaled kinetic spectrum Ekin(k)/R*\n"
        "Gamma : float, optional (default: 4.0/3.0)\n    Mean adiabatic index"
        " Gamma = w_bar/e_bar of the fluid in the stable phase\n"
        "\nReturns\n-------\n"
        "OmegaK : float\n    kinetic energy density parameter OmegaK"
    );

    m.def(
        "kinetic_energy_density_parameter",
        [](
            double xi_wall,
            py::array_t<double> chi_arr, 
            py::array_t<double> A2_arr, 
            double Gamma
        ) { 
            return kinetic_energy_density_parameter(
                xi_wall, np2vec(chi_arr), np2vec(A2_arr), Gamma
            );
        },
        py::arg("xi_wall"), 
        py::arg("chi_arr"), 
        py::arg("A2_arr"), 
        py::arg("Gamma")=4.0/3.0,
        "Total kinetic energy density parameter OmegaK.\n\n"
        "This is the total kinetic energy density OmegaK = rho_kin/rho_crit as"
        " a fraction of the critical energy density.\n"
        "The integral over T~ in Eq. (4.33) has been carried out analytically,"
        " assuming that the nucleation history is either exponential or"
        " simultaneous, as in the step from Eq. (4.33) to (4.34) in"
        " [Hindmarsh+ (2019)].\n"
        "\nParameters\n----------\n"
        "chi_arr : numpy.ndarray, shape (N,)\n"
        "    Values of chi = T q = T~ q/beta at which |A|^2 is provided\n"
        "A2_arr : numpy.ndarray, shape (N,)\n    Precomputed values of the"
        " shape function |A|^2 used in the integration\n"
        "xi_wall : float\n    Wall velocity xi_w\n"
        "\nReturns\n-------\n"
        "OmegaK : float\n    kinetic energy density parameter OmegaK"
    );

    m.def(
        "anisotropic_stress_power_spectrum_amplitude",
        [](py::array_t<double> kR, py::array_t<double> zeta_kin) { 
            return anisotropic_stress_power_spectrum_amplitude(
                np2vec(kR), np2vec(zeta_kin)
            ); 
        },
        py::arg("kR"), py::arg("zeta_kin"),
        "Normalized amplitude script_C of the anisotropic stress power spectrum"
        " k E_Pi.\n\nThis is the constant script_C defined in Eq. (46) of"
        " [Roper+ (2024)]. It describes the normalized amplitude of the"
        " anisotropic stress power spectrum [cf. Eq. (44) of [Roper+ (2024)]]\n"
        "    k E_Pi(tau_1,tau_2,k) = 2 w^2 K^3 (Omega_K/script_K)^2"
        " script_C zeta_Pi(tau-, K)\n"
        "where K = k R*.\n"
        "\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n    Momentum k times bubble separation"
        " R* at which the kinetic spectral shape is evaluated\n"
        "zeta_kin : numpy.ndarray, shape (N,)\n"
        "    Kinetic spectral shape zeta_kin(k R*)\n"
        "\nReturns\n-------\n"
        "script_C : float\n    Normalized amplitude script_C of the anisotropic"
        " stress power spectrum k E_Pi"
    );

    m.def(
        "anisotropic_stress_correlator_integration",
        [](
            py::array_t<double> kR, 
            utils::UniformLinearInterpolator& E_kin_over_R, 
            AnisotropicStressCorrelatorSettings settings
        ) {
            return vec2np(anisotropic_stress_correlator_integration(
                np2vec(kR), E_kin_over_R, settings
            ));
        },
        py::arg("kR"), py::arg("E_kin_over_R"), py::arg("settings"),
        "Calculate the nested momentum integral of the anisotropic stress UETC"
        " E_Pi or gravitational wave spectrum OmegaGW*.\n\n"
        "This function calculates the nested momentum (and time) integrals for"
        " different choices of the time dependence function Delta. The choice"
        " of Delta is controlled by the universe and freq_approx parameters of"
        " the AnisotropicStressCorrelatorSettings settings, and the prefactor"
        " w_bar is set from settings.c_w.\n"
        "If settings.universe is set to ExpansionHistory.NONE (settings.c_w to"
        " 2 w_bar where w_bar is the mean enthalpy density), the function"
        " returns the anisotropic stress (auto-correlation) power spectrum"
        " k E_Pi(k) in Eq. (44) of [Roper+ (2024)].\n"
        "For ExpansionHistory.RADIATION_DOMINATION or"
        " ExpansionHistory.STATIC_UNIVERSE (and settings.c_w =  3 Gamma^2), the"
        " function returns the gravitational wave spectrum OmegaGW(k) at"
        " production, i.e. Eq. (51) of [Roper+ (2024)] setting T_GW=1, with"
        " Delta assuming an radiation dominated or static universe,"
        " respectively, and using the approximation specified in"
        " settings.freq_approx.\n"
        "\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n    Momentum k times bubble separation"
        " R* at which the anisotropic stress UETC is evaluated\n"
        "E_kin_over_R : UniformLinearInterpolator\n    Interpolating function"
        " for the scaled kinetic spectrum Ekin/R* as a function of log(p R_*)\n"
        "settings : AnisotropicStressCorrelatorSettings\n"
        "    AnisotropicStressCorrelatorSettings object containing the settings"
        " for the momentum and time integration\n"
        "\nReturns\n-------\n"
        "K3_zeta : numpy.ndarray, shape (N,)\n"
        "    Integrated anisotropic stress UETC or GW spectrum"
    );

    m.def(
        "anisotropic_stress_power_spectrum",
        [](
            py::array_t<double> kR, 
            utils::UniformLinearInterpolator& E_kin_over_R, 
            double w_bar, 
            size_t N_p, 
            size_t N_z
        ) {
            return vec2np(anisotropic_stress_power_spectrum(
                np2vec(kR), E_kin_over_R, w_bar, N_p, N_z
            ));
        },
        py::arg("kR"), 
        py::arg("E_kin_over_R"), 
        py::arg("w_bar")=1.0, 
        py::arg("N_p")=5000, 
        py::arg("N_z")=100,
        "Anisotropic stress auto-correlation power spectrum k E_Pi.\n\n"
        "This is the power spectrum of the anisotropic stress auto-correlation,"
        " k E_Pi(k) = k E_Pi(tau,tau,k), Eq. (47) of [Roper+ (2024)] with"
        " tau-=0, i.e. the integration in anisotropic_stress_correlator_"
        "integration with Delta = 1 (settings.universe = ExpansionHistory.NONE)"
        " and c_w = 2 w_bar^2 (settings.c_w = 2.0).\n"
        "\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n    Momentum k times bubble separation"
        " R* at which the anisotropic stress UETC is evaluated\n"
        "E_kin_over_R : UniformLinearInterpolator\n    Interpolating function"
        " for the scaled kinetic spectrum Ekin/R* as a function of log(p R_*)\n"
        "w_bar : float, optional (default: 1.0)\n"
        "    Mean enthalpy density w_bar (far ahead of the wall)\n"
        "N_p : integer, optional (default: 5000)\n"
        "    Number of points in the integration over p R*\n"
        "N_z : integer, optional (default: 100)\n"
        "    Number of points in the integration over z\n"
        "\nReturns\n-------\n"
        "k_E_Pi : numpy.ndarray, shape (N,)\n    Anisotropic stress"
        " auto-correlator power spectrum k E_Pi(tau,tau,k)"
    );

    m.def(
        "anisotropic_stress_power_spectral_shape",
        [](
            py::array_t<double> k_E_Pi, 
            double script_C, 
            double E_kin_max_over_R, 
            double w_bar
        ) {
            return vec2np(anisotropic_stress_power_spectral_shape(
                np2vec(k_E_Pi), script_C, E_kin_max_over_R, w_bar
            ));
        },
        py::arg("k_E_Pi"), 
        py::arg("script_C"), 
        py::arg("E_kin_max_over_R"), 
        py::arg("w_bar")=1.0,
        "Spectral shape of the anisotropic stress UETC power spectrum (k R*)^3"
        " zeta_Pi(k R*).\n\n"
        "This is the combination K^3 zeta_Pi(K) with script_C zeta_Pi(K)"
        " defined in Eq. (47) of [Roper+ (2024)]. It describes the spectral"
        " shape of the power spectrum of the anisotropic stress UETC k E_Pi."
        " The constant script_C can be calculated with"
        " anisotropic_stress_power_spectrum_amplitude(qR, zeta_kin).\n"
        "\nParameters\n----------\n"
        "k_E_Pi : numpy.ndarray, shape (N,)\n"
        "    Anisotropic stress auto-correlator k E_Pi\n" 
        "script_C : float\n    Amplitude parameter script_C of the anisotropic"
        " stress power spectrum\n"
        "E_kin_max_over_R : float\n    Maximal value Ekin* of the kinetic"
        " spectrum Ekin(k) normalized to the average bubble size R*\n"
        "w_bar : float, optional (default: 1.0)\n"
        "    Mean enthalpy density w_bar (far ahead of the wall)\n"
        "\nReturns\n-------\n"
        "K3_zeta_Pi : numpy.ndarray, shape (N,)\n    Spectral shape of the"
        " anisotropic stress UETC power spectrum (k R*)^3 zeta_Pi(k R*)"
    );  

    m.def(
        "infinite_duration_gw_power_spectrum",
        [](
            py::array_t<double> kR, 
            utils::UniformLinearInterpolator& E_kin_over_R, 
            double Upsilon_H_R, 
            size_t N_p, 
            double c_w
        ) {
            return vec2np(infinite_duration_gw_power_spectrum(
                np2vec(kR), E_kin_over_R, Upsilon_H_R, N_p, c_w
            ));
        },
        py::arg("kR"), 
        py::arg("E_kin_over_R"), 
        py::arg("Upsilon_H_R"), 
        py::arg("N_p")=5000, 
        py::arg("c_w")=16.0/3.0,
        "Approximated gravitational wave power spectrum OmegaGW* in the"
        " infinite duration limit.\n\n"
        "This is the gravitational wave power spectrum OmegaGW* at production"
        " [cf. anisotropic_stress_correlator_integration] in the inifinite"
        " source duration limit, Eq. (B3) of [Roper+ (2024)] (with T_GW=1)."
        " In the k >> 1/delta_tau_fin limit, the time integration factor Delta"
        " becomes a Dirac-delta distribution, simplifying the momentum"
        " integration. The suppression factor Upsilon due to the finite"
        " lifetime of the soundwave source is given by Upsilon = H*"
        " delta_tau_fin/(1 + H* delta_tau_fin) in a radiation-dominated"
        " universe and H* delta_tau_fin in a static universe. This corresponds"
        " to P_GW in Eq. (3.48) of [Hindmarsh+ (2019)].\n"
        "\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n    Momentum k times bubble separation"
        " R* at which the spectrum is evaluated\n"
        "E_kin_over_R : UniformLinearInterpolator\n    Interpolating function"
        " for the scaled kinetic spectrum Ekin/R* as a function of log(p R_*)\n"
        "Upsilon_H_R : float\n"
        "    Finite source lifetime suppression factor, Upsilon, times the mean"
        " bubble separation normalized to the Hubble radius, H* R*\n"
        "N_p : integer, optional (default: 5000)\n"
        "    Number of points in the integration over p R*\n"
        "c_w : float, optional (default: 16.0/3.0)\n"
        "    Prefactor c_w in anisotropic_stress_correlator_integration, should"
        " be set to c_w = 3 Gamma^2\n"
        "\nReturns\n-------\n"
        "OmegaGW : numpy.ndarray, shape (N,)\n    GW power spectrum OmegaGW* at"
        " production in the k delta_tau_fin >> 1 limit"
    );

    m.def(
        "gravitational_wave_power_spectrum",
        [](
            py::array_t<double> kR, 
            utils::UniformLinearInterpolator& E_kin_over_R,
            double H_delta_tau,
            double H_R,
            ExpansionHistory expansion,
            MomentumIntegration integration,
            double Gamma,
            size_t N_p,
            size_t N_z
        ) { return vec2np(gravitational_wave_power_spectrum(
            np2vec(kR), 
            E_kin_over_R,
            H_delta_tau, 
            H_R, 
            expansion, 
            integration, 
            Gamma, 
            N_p, 
            N_z
        )); },
        py::arg("kR"), 
        py::arg("E_kin_over_R"), 
        py::arg("H_delta_tau"), 
        py::arg("H_R"), 
        py::arg("expansion")=ExpansionHistory::RADIATION_DOMINATION,
        py::arg("integration")=MomentumIntegration::FULL, 
        py::arg("Gamma")=4.0/3.0, 
        py::arg("N_p")=5000, 
        py::arg("N_z")=100,
        "Gravitational wave power spectrum OmegaGW*(k) at production.\n\n"
        "This is Eq. (51) of [Roper+ (2024)] (with T_GW=1).\n"
        "\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n    Momentum k times bubble separation"
        " R* at which the spectrum is evaluated\n"
        "E_kin_over_R : UniformLinearInterpolator\n    Interpolating function"
        " for the scaled kinetic spectrum Ekin/R* as a function of log(p R_*)\n"
        "H_delta_tau : float\n"
        "    Source duration in units of the Hubble time, H* delta_tau_fin\n"
        "H_R : float\n"
        "    Mean bubble separation in units of the Hubble size, H* R*\n"
        "expansion : ExpansionHistory, optional (default:"
        " ExpansionHistory.RADIATION_DOMINATION\n    Assumed ExpansionHistory\n"
        "integration : MomentumIntegration, optional (default:"
        " MomentumIntegration.FULL)\n    Method for calculating the momentum"
        " integrals in anisotropic_stress_correlator_integration\n"
        "Gamma : float, optional (default: 4/3)\n    Mean adiabatic index"
        " Gamma = w/e of the fluid in the stable phase\n"
        "N_p : integer, optional (default: 5000)\n"
        "    Number of points in the integration over p R*\n"
        "N_z : integer, optional (default: 100)\n"
        "    Number of points in the integration over z\n"
        "\nReturns\n-------\n"
        "OmegaGW : numpy.ndarray, shape (N,)\n"
        "    GW power spectrum OmegaGW*(k) at production"
    );

    m.def(
        "Delta_rad",
        &Delta_rad,
        py::arg("H_delta_tau"), 
        py::arg("k_tau"), 
        py::arg("p_tau"), 
        py::arg("pt_tau"),
        "Time integral function Delta(delta_tau_fin, k, p, p~) for a radiation"
        " dominated universe.\n\n"
        "This is the function Delta that contains the integration over"
        " conformal times tau1 and tau2 of the time dependence in the Green's"
        " function and the stationary UETC (cf. Eq. (52) of [Roper+ (2024)]),"
        " assuming a radiation dominated universe.\n\n"
        "We use the dimensionless variables q tau* = q/H* for the wave numbers"
        " q=k,p,p~ and H delta_tau_fin = delta_tau_fin/tau* as inputs.\n"
        "Note that, since we use a* = a(tau*) = 1, we have H* = 1/tau*, so that"
        " q tau* = (q R*)/(H* R*).\n"
        "\nParameters\n----------\n"
        "H_delta_tau : float\n    Duration of the gravitational wave source"
        " normalized to the Hubble time H* delta_tau_fin\n"
        "k_tau : float\n    Dimensionless wave number k tau*\n"
        "p_tau : float\n    Dimensionless wave number p tau*\n"
        "pt_tau : float\n    Dimensionless wave number p~ tau*\n"
        "\nReturns\n-------\n"
        "Delta : float\n"
        "    Value of the time integral Delta(delta_tau_fin, k, p, p~)"
    );

    m.def(
        "Delta_mn",
        &Delta_mn,
        py::arg("H_delta_tau"), py::arg("pmn_tau"),
        "Intermediate function Delta_mn(delta_tau_fin, p_mn) for the time"
        " integral in a radiation-dominated universe.\n\n"
        "This is the intermediate function Delta_mn(delta_tau_fin, p_mn) used"
        " in the calculation of the time integral Delta_rad, cf. Eq. (56) of"
        " [Roper+ (2024)].\n\n"
        "We use the dimensionless variables pmn tau* = pmn/H* and H"
        " delta_tau_fin = delta_tau_fin/tau* as inputs.\n"
        "\nParameters\n----------\n"
        "H_delta_tau : float\n    Duration of the gravitational wave source"
        " normalized to the Hubble time H* delta_tau_fin\n"
        "pmn_tau : float\n    Dimensionless wave number p_mn tau*\n"
        "\nReturns\n-------\n"
        "Delta_mn : float\n    Value of the intermediate time integral"
        " contribution Delta_mn(delta_tau_fin, p_mn)"
    );

    m.def(
        "Delta_static",
        &Delta_static,
        py::arg("H_delta_tau"),
        py::arg("k_tau"), 
        py::arg("p_tau"), 
        py::arg("pt_tau"),
        "Time integral function Delta(delta_tau_fin, k, p, p~) for a static"
        " universe.\n\n"
        "This is the function Delta that contains the integration over"
        " conformal times tau1 and tau2 of the time dependence in the Green's"
        " function and the stationary UETC (cf. Eq. (52) of [Roper+ (2024)]),"
        " assuming a static universe.\n\n"
        "We use the dimensionless variables q tau* = q/H* for the wave numbers"
        " q=k,p,p~ and H delta_tau_fin = delta_tau_fin/tau* as inputs.\n"
        "Note that, since we use a* = a(tau*) = 1, we have H* = 1/tau*, so that"
        " q tau* = (q R*)/(H* R*).\n"
        "\nParameters\n----------\n"
        "H_delta_tau : float\n    Duration of the gravitational wave source"
        " normalized to the Hubble time H* delta_tau_fin\n"
        "k_tau : float\n    Dimensionless wave number k tau*\n"
        "p_tau : float\n    Dimensionless wave number p tau*\n"
        "pt_tau : float\n    Dimensionless wave number p~ tau*\n"
        "\nReturns\n-------\n"
        "Delta : float\n    Value of the time integral"
        " Delta_flat(delta_tau_fin, k, p, p~) for a static universe"
    );

    m.def(
        "Delta_mn_static",
        &Delta_mn_static,
        py::arg("H_delta_tau"), py::arg("pmn_tau"),
        "Intermediate function Delta_mn_flat(delta_tau_fin, p_mn) for the time"
        " integral in a static universe.\n\n"
        "This is Eq. (59) of [Roper+ (2024)].\n\n"
        "We use the dimensionless variables pmn tau* = pmn/H* and H"
        " delta_tau_fin = delta_tau_fin/tau* as inputs.\n"
        "\nParameters\n----------\n"
        "H_delta_tau : float\n    Duration of the gravitational wave source"
        " normalized to the Hubble time H* delta_tau_fin\n"
        "pmn_tau : float\n    Dimensionless wave number p_mn tau*\n"
        "\nReturns\n-------\n"
        "Delta_mn : float\n"
        "    Value of the time integral contribution Delta_mn_flat"
        "(delta_tau_fin, p_mn) neglecting the expansion of the Universe"
    );

    m.def(
        "Delta_low_freq_rad",
        &Delta_low_freq_rad,
        py::arg("H_delta_tau"), py::arg("p_tau"),
        "Time integral function in the low frequency limit"
        " Delta0(delta_tau_fin, p) = lim(k->0) Delta(delta_tau_fin,k,p,p~).\n\n"
        "This is the low-frequency limit Delta0 of the integration over"
        " conformal times tau1 and tau2 of the time dependence in the Green's"
        " function and the stationary UETC, cf. Eq. (61) of [Roper+ (2024)],"
        " assuming radiation domination.\n\n"
        "We use the dimensionless variables p tau* = p/H* and H delta_tau_fin"
        " = delta_tau_fin/tau* as inputs.\n"
        "Note that, since we use a* = a(tau*) = 1, we have H* = 1/tau*, so that"
        " p tau* = (p R*)/(H* R*).\n"
        "\nParameters\n----------\n"
        "H_delta_tau : float\n    Duration of the gravitational wave source"
        " normalized to the Hubble time H* delta_tau_fin\n"
        "p_tau : float\n    Dimensionless wave number p tau*\n"
        "\nReturns\n-------\n"
        "Delta0 : float\n"
        "    Value of the time integral Delta0 in the low frequency limit"
    );

    m.def(
        "Delta_low_freq_static",
        &Delta_low_freq_static,
        py::arg("H_delta_tau"), py::arg("p_tau"),
        "Time integral function in the low frequency limit in a static"
        " universe, Delta0_flat(delta_tau_fin, p) = lim(k->0)"
        " Delta_flat(delta_tau_fin,k,p,p~).\n\n"
        "This is the low-frequency limit Delta0_flat of the integration over"
        " conformal times tau1 and tau2 of the time dependence in the Green's"
        " function and the stationary UETC, cf. Eq. (61) of [Roper+ (2024)],"
        " assuming a static universe.\n\n"
        "We use the dimensionless variables p tau* = p/H* and H delta_tau_fin"
        " = delta_tau_fin/tau* as inputs.\n"
        "Note that, since we use a* = a(tau*) = 1, we have H* = 1/tau*, so that"
        " p tau* = (p R*)/(H* R*).\n"
        "\nParameters\n----------\n"
        "H_delta_tau : float\n    Duration of the gravitational wave source"
        " normalized to the Hubble time H* delta_tau_fin\n"
        "p_tau : float\n    Dimensionless wave number p tau*\n"
        "\nReturns\n-------\n"
        "Delta0 : float\n    Value of the time integral Delta0_flat in the low"
        " frequency limit in a static universe"
    );

    m.def(
        "weighted_average_Delta_0",
        [](
            py::array_t<double> kR, 
            py::array_t<double> E_kin_over_R,  
            double H_delta_tau, 
            double H_R, 
            ExpansionHistory expansion
        ) { return weighted_average_Delta_0(
            np2vec(kR), np2vec(E_kin_over_R), H_delta_tau, H_R, expansion
        ); },
        py::arg("kR"), 
        py::arg("E_kin_over_R"), 
        py::arg("H_delta_tau"), 
        py::arg("H_R"), 
        py::arg("expansion")=ExpansionHistory::RADIATION_DOMINATION,
        "Weighted average Delta0~ of Delta0 weighted with the kinetic spectrum"
        " Ekin.\n\n"
        "This is the weighted average of the low-frequency integral function"
        " Delta0 with the kinetic spectrum Ekin, parametrizing the dependence"
        " of the GW amplitude on the source duration delta_tau_fin, cf. Eq."
        " (63) of [Roper+ (2024)].\n"
        "\nParameters\n----------\n"
        "kR : numpy.ndarray, shape (N,)\n"
        "    Momentum k times bubble separation R*\n"
        "E_kin_over_R : numpy.ndarray, shape (N,)\n"
        "    Scaled kinetic spectrum Ekin(k)/R*\n"
        "H_delta_tau : float\n    Duration of the gravitational wave source"
        " normalized to the Hubble time H* delta_tau_fin\n"
        "H_R : float\n"
        "    Mean bubble separation in units of the Hubble size, H* R*\n"
        "expansion : ExpansionHistory, optional (default:"
        " ExpansionHistory.RADIATION_DOMINATION\n    Assumed ExpansionHistory\n"
        "\nReturns\n-------\n"
        "Delta0_tilde : float:\n    Weighted average source lifetime dependence"
        " Delta0~(delta_tau_fin, R*)"
    );

    m.def(
        "normalized_GW_spectrum",
        [](
            py::array_t<double> Omega_GW, 
            double script_C, 
            double E_kin_max_over_R, 
            double Delta_tilde_0, 
            double Gamma
        ) {
            return vec2np(normalized_GW_spectrum(
                np2vec(Omega_GW), 
                script_C, 
                E_kin_max_over_R, 
                Delta_tilde_0, 
                Gamma
            ));
        },
        py::arg("Omega_GW"), 
        py::arg("script_C"), 
        py::arg("E_kin_max_over_R"), 
        py::arg("Delta_tilde_0"), 
        py::arg("Gamma")=4.0/3.0,
        "Normalized gravitational wave spectrum (k R*)^3 zeta_GW(k R*).\n\n"
        "This is the normalized gravitational wave spectrum in Eq. (93) of"
        " [Roper+ (2024)].\n"
        "\nParameters\n----------\n"
        "Omega_GW : numpy.ndarray, shape (N,)\n"
        "    Gravitational wave power spectrum OmegaGW*(k) at production,"
        " cf. gravitational_wave_power_spectrum\n"
        "script_C : float\n"
        "    Amplitude parameter script_C of the anisotropic stress power"
        " spectrum, cf. anisotropic_stress_power_spectrum_amplitude\n"
        "E_kin_max_over_R : float\n    Maximal value Ekin* of the kinetic"
        " spectrum Ekin normalized to the average bubble size R*\n"
        "Delta_tilde_0 : float\n    Weighted average source lifetime dependence"
        " Delta0~, cf. weighted_average_Delta_0\n"
        "Gamma : float, optional (default: 4/3)\n    Mean adiabatic index"
        " Gamma = w/e of the fluid in the stable phase\n"
        "\nReturns\n-------\n"
        "K3_zeta_GW : numpy.ndarray, shape (N,)\n"
        "    Normalized gravitational wave spectrum K^3 zeta_GW"
    );

    m.def(
        "spectral_modification_function",
        [](py::array_t<double> zeta_GW, py::array_t<double> zeta_Pi) { 
            return vec2np(spectral_modification_function(
                np2vec(zeta_GW), np2vec(zeta_Pi)
            )); 
        },
        py::arg("zeta_GW"),  py::arg("zeta_Pi"),
        "Spectral modification function Delta~ = zeta_GW/zeta_Pi.\n\n"
        "This functions parametrized the modification of the spectral shape of"
        " the graviational wave spectrum due to the finite source duration,"
        " Eq. (95) of [Roper+ (2024)].\n"
        "\nParameters\n----------\n"
        "zeta_GW : numpy.ndarray, shape (N,)\n    Normalized gravitational wave"
        " spectrum zeta_GW (or K^3 zeta_GW) [cf. normalized_GW_spectrum]\n"
        "zeta_Pi : numpy.ndarray, shape (N,)\n"
        "    Normalized anisotropic stress spectrum zeta_Pi (or K^3 zeta_Pi)"
        " [cf. anisotropic_stress_power_spectral_shape]\n"
        "\nReturns\n-------\n"
        "Delta_tilde : numpy.ndarray, shape (N,)\n"
        "    Spectral modification function Delta~ = zeta_GW/zeta_Pi"
    );

}