# Inverse PT Calculator

[![license](https://img.shields.io/github/license/eric-madge/inverse_pt)](
  https://github.com/eric-madge/inverse_pt/blob/main/LICENSE
)
[![arXiv](https://img.shields.io/badge/arXiv-2510.21439-%23B31B1B)](
  https://arxiv.org/abs/2510.21439
)

A C++ code for the calculation fluid profiles of direct and inverse phase 
transitions and for computing the corresponding gravitational wave spectrum 
within the sound shell model.

Contact: [eric.madgepimentel@uam.es](mailto:eric.madgepimentel@uam.es)

## Table of Contents

- [Description](#description)
- [Attribution](#attribution)
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)
- [References](#references)

## Description

This code provides the C++ library `libinverse_pt_calculator`. It is structured
into four main header files:

1. The header `include/inverse_pt/profile_calculator.hpp` defines the namespace
  `inverse_pt::profile_calculator`, which provides the functions to calculate
  the velocity, enthalpy, energy and pressure profiles of the fluid for given
  values of the wall velocity $\xi_w$ and transition strength $\alpha_N$. 
2. The header `include/inverse_pt/fluid_profile.hpp` provides the class
  `inverse_pt::FluidProfile`. This class acts as a simple interface to the 
  functions in `inverse_pt::profile_calculator`, automatically passing the 
  results from one routine to the next.
3. The header `include/inverse_pt/sound_shell_model.hpp` defines the namespace
  `inverse_pt::sound_shell_model` with the functions to compute the
  gravitational wave spectrum and various related quantities within the sound
  shell model.
4. The header `include/inverse_pt/sound_shell_spectrum.hpp` provides the class
  `inverse_pt::SoundShellSpectrum`. This is an object-oriented interface to the
  functions in `inverse_pt::sound_shell_model`.

The `include` folder furhermore contains the following files:
5. The header `include/inverse_pt/utils.hpp` and the corresponding namespace
  `inverse_pt::utils` contain helper functions used inside the code.
6. The header `include/inverse_pt/settings.hpp` defines the namespace
  `inverse_pt::settings` which includes global settings, such as a flag `debug`
  for additional debug prints, or a `check_mode` switch to control handling of 
  non-critical consistency checks.
7. The header `include/inverse_pt/constants.hpp` defines the namespace
  `inverse_pt::constants` containing the physical and  mathematical constants 
  used in the code.

The Python module `inverse_pt_calculator` provides the submodules 
`inverse_pt_calculator.profile_calculator` and 
`inverse_pt_calculator.sound_shell_model` as well as the classes
`inverse_pt_calculator.FluidProfile` and 
`inverse_pt_calculator.SoundShellSpectum` which interface to the corresponding
C++ versions.

## Attribution

If you use this code in your work, please cite
<a href="#Barni+ (2025b)">Barni, Blasi, Madge and Vanvlasselaer (2025)
[Barni+ (2025b)]</a>.  
```bibtex
@article{Barni:2025gnm,
    author = "Barni, Giulio and Blasi, Simone and Madge, Eric and Vanvlasselaer, Miguel",
    title = "{Gravitational waves from the sound shell model: direct and inverse phase transitions in the early Universe}",
    eprint = "2510.21439",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "DESY-25-144, IFT-UAM/CSIC-25-118",
    month = "10",
    year = "2025"
}
```

## Installation

### Prerequisites 

To install and run the code, you need the following prerequisites.

#### Build Tools

- CMake ≥ 3.24
- C++17 compatible compiler (tested with Clang/GCC)
- Python ≥ 3.6 (for building the Python interface if `BUILD_PYTHONLIB=ON`)

#### Dependencies

- [Boost](https://www.boost.org/)
- [GSL (GNU Scientific Library)](https://www.gnu.org/software/gsl/)
- [OpenMP](https://www.openmp.org/)
- [pybind11](https://github.com/pybind/pybind11) 
  (for building the Python interface if `BUILD_PYTHONLIB=ON`)
- [Doxygen](https://www.doxygen.nl/) 
  (for building the documentation if `BUILD_DOCS=ON`)

### Basic Instructions

The `inverse_pt_calculator` code can be installed using cmake with the
following commands:
```bash
cmake -S<root_dir> -B<build_folder> -DCMAKE_INSTALL_PREFIX=<install_path>
cmake --build <build_folder> -j <N>
cmake --install <build_folder>
```
or
```bash
mkdir <build_folder> && cd <build_folder>
cmake <root_dir> -DCMAKE_INSTALL_PREFIX=<install_path>
make -j <N>
make install
```
where `<root_dir>` is the root directory of the code (i.e. the one containing
the `CMakeLists.txt` file), `<build_folder>` is the folder in which the project
will be built,`<install_path>` is the path into which the code will be 
installed, and `<N>` is the number of cores when building in parallel.

For instance, if the project root directory is `inversePT`, you want to build
the code in the folder `inversePT/build` on 8 cores and install it into
`inversePT/install`, you can go to `inversePT` and run
```bash
cmake -S. -Bbuild -DCMAKE_INSTALL_PREFIX=install
cmake --build build -j 8
cmake --install build
```

After installation, the `<install_path>` will contain the following files
- `lib/libinverse_pt_calculator.so` or `lib/libinverse_pt_calculator.dylib`  
  The shared `inverse_pt_calculator` library to which you can link your code.
- `include/inverse_pt/*.hpp`  
  Header files for the C++ library that you need to include to use the library.
- `bin/test` (if `BUILD_TEST=ON`)  
  A test executable that calculates the fluid profiles and gravitational wave
  spectrum for given input values.
- `lib/python<python version>/site-packages/inverse_pt_calculator.cpython-<python version>-<platform>.so`
  (if `BUILD_PYTHONLIB=ON`)  
  Python interface to the library.
- `share/inverse_pt_calculator/docs/` (if `BUILD_DOCS=ON`)    
  HTML and LaTeX documentation generated with Doxygen.

### Dependencies

#### Compiler

This project requires a C++17-compatible compiler, such as:
- GCC >= 7
- Clang >= 6
- Apple Clang >= 10

You can select a specific compiler using environment variables:
```bash
export CXX=/path/to/g++
```
or pass it directly to CMake:
```bash
cmake -DCMAKE_CXX_COMPILER=/path/to/clang++ ..
```

#### OpenMP

OpenMP is used for parallel computation. CMake will automatically find it if 
available. If not detected:
- Install the relevant OpenMP development package (e.g. `libomp-dev` on 
Debian/Ubuntu, or `libomp` via MacPorts or Homebrew on macOS).
- If you have a compiler with working OpenMP but CMake complains about not 
  finding OpenMP, make sure that it is using the correct compiler (see 
  [above](#compiler)). E.g., to use the compiler that is invoked when you 
  call `g++` in the terminal, you can use the option
  `-DCMAKE_CXX_COMPILER=$(which g++)`.
- If necessary, set: `-DOpenMP_CXX_FLAGS="-fopenmp"`,
  `-DOpenMP_CXX_LIB_NAMES="omp"` and `-DOpenMP_omp_LIBRARY=/path/to/libomp.a`


#### GSL

The GNU Scientific Library is required for evaluating the cosine and sine
integrals in the integration of the time dependence of the anisotropic stress
unequal time correlator. If it is not detected automatically, you can help
CMake find it by setting `-DGSL_ROOT_DIR=/path/to/gsl`. If GSL is installed via 
a package manager, you can use `pkg-config gsl --variable=prefix` to locate it.

#### Boost

The `inverse_pt_calculator` code uses the Boost library for root finding. If 
Boost is not found automatically, you can, depending on your CMake version, 
provide one of the following hints:
- For CMake >= 3.30 (i.e. with policy [`CMP0167`](
    https://cmake.org/cmake/help/latest/policy/CMP0167.html#policy:CMP0167
  ) set to `NEW`), you can provide the path to Boost's `BoostConfig.cmake` CMake
  configuration file via the `-DBoost_DIR=/path/to/boost/cmake` option, e.g. 
  for Boost-1.81 installed with MacPorts:
  ```bash
  cmake -DBoost_DIR=/opt/local/libexec/boost/1.81/lib/cmake/Boost-1.81.0
  ```
- For older CMake versions (i.e. with policy [`CMP0167`](
    https://cmake.org/cmake/help/latest/policy/CMP0167.html#policy:CMP0167
  ) set to `OLD`), you can provide the path to Boost's installation prefix via
  the `-DBoost_ROOT=/path/to/boost/root` option, e.g.
  ```bash
  cmake -DBoost_ROOT=/opt/local/libexec/boost/1.81/
  ```
  You might also have to provide `BOOST_INCLUDEDIR` and `BOOST_LIBRARYDIR`.

#### Python3

Python3 is an optional dependence required for building the Python interface.
If CMake does not find it automatically, or if you want to uses a specific 
Python version installed on your system, you can provide the path to the 
`-DPython3_ROOT_DIR=/path/to/python/root`. You can obtain the path from the
corresponding `python3-config` executable by
```bash
python3-config --prefix
```

#### pybind11

Pybind11 is optional but required if the Python interface is built. You can
provide the paths to the packages CMake configuration file with
`-Dpybind11_DIR=/path/to/pybind11/cmake`, which can be found via
```bash
python3 -m pybind11 --cmakedir
```
If CMake still complains about not finding pybind11, make sure that it is using
the corresponding Python3 interpreter (cf. [Python3](#python3)).

#### Doxygen

Doxygen is optional and required only for generating documentation. This step
can be skipped by setting `-DBUILD_DOCS=OFF`. If your doxygen executable is not
in a standard path, you can set the path explicitly with
`-DDoxygen::doxygen=/path/to/doxygen`.

### Configuration

Pass `-D<OPTION>=<ON|OFF>` to the cmake configureation command
(`cmake -S<root_dir> -B<build_folder>` or `cmake <root_dir>`) to switch the
option `<OPTION>` on or off.

The following options control which components of the code are built:

- `BUILD_PYTHON_LIB`
  Whether to build (and install) the Python interface. Building the interface
  requires Python and pybind11 (which in turn requires Cython). This option is
   set to `ON` by default.
- `BUILD_DOCS`
  Whether to build (and install) the documentation. Building the documentation
  requires Doxygen. This option is set to `ON` by default.
- `BUILD_TEST`
  Whether to build (and install) the test executable. This option is set to
  `ON` by default.

Furthermore, you can further control the building process using
- `CMAKE_POSITION_INDEPENDENT_CODE`
  Generate position independent code (`-fPIC` flag). This option is set to `ON`
  by default.
- `CMAKE_INTERPROCEDURAL_OPTIMIZATION`
  Use link-time optimization (`-flto` flag). This option is set to `ON` by 
  default.
- `OPTIMIZE_NATIVE`
  Optimize the code for the specific machine architecture (`-march=native`
  flag). This option is set to `ON` by default.
- `DEVELOPMENT_WARNINGS`
  Show compilation warnings (`-Wall` flag, for development). This option is set
  to `OFF` by default. Use the `--verbose` option in the build process to see 
  the warnings.
- `VECTORIZATION_REPORTS`
  Show vectorization status reports (`--Rpass=loop-vectorize` and
  `-Rpass-missed=loop-vectorize` flags, for development). This option is set to
  `OFF` by default. Use the `--verbose` option in the build process to see the
  reports.
- `QUIET_DOXYGEN`
  Ony display warnings and errors when creating the documentation. This option
  is set to `ON` by default.

## Usage

No usage information yet.

Please consult the html documentation in the `docs` folder for the C++ library
or use the `help` function in Python. 

You can import the Python module by
```
import sys
sys.path.append('<install_path>/lib/python<python version>/site-packages/')
import inverse_pt_calculator as iptc
```
and then `help(iptc)` for further information on the package content.

## License

The inverse PT calculator code is licensed under the GNU General Public License,
version 3 or later. See the accompanying [LICENSE file](LICENSE) for further
details.

## References

<a id="Espinosa+ (2010)">[Espinosa+ (2010)]</a>    
  J. R. Espinosa and J. M. No,  
  *Energy Budget of Cosmological First-order Phase Transitions*,  
  [JCAP **06** (2010), 028](https://doi.org/10.1088/1475-7516/2010/06/028) 
  [[arXiv:1004.4187 [hep-ph]](https://arxiv.org/abs/1004.4187)].

<a id="Hindmarsh (2018)">[Hindmarsh (2018)]</a>    
  M.Hindmarsh,  
  *Sound shell model for acoustic gravitational wave production at a first-order phase transition in the early Universe*,  
  [Phys. Rev. Lett. **120** (2018) no. 7, 071301](
    https://doi.org/10.1103/PhysRevLett.120.071301
  ) 
  [[arXiv: 1608.04735 [astro-ph.CO]](https://arxiv.org/abs/1608.04735)].

<a id="Hindmarsh+ (2019)">[Hindmarsh+ (2019)]</a>    
  M. Hindmarsh and M. Hijazi  
  *Gravitational waves from first order cosmological phase transitions in the Sound Shell Model*,  
  [JCAP **12** (2019), 062](https://doi.org/10.1088/1475-7516/2019/12/062) 
  [[arXiv:1909.10040 [astro-ph.CO]](https://arxiv.org/abs/1909.10040)].

<a id="Roper+ (2024)">[Roper+ (2024)]</a>    
  A. Roper Pol, S. Procacci and C. Caprini,  
  *Characterization of the gravitational wave spectrum from sound waves within
  the sound shell model*,  
  [Phys. Rev. D **109** (2024) 6, 063531](
    https://doi.org/10.1103/PhysRevD.109.063531
  ) 
  [[arXiv:2308.12943 [gr-qc]](https://arxiv.org/abs/2308.12943)].

<a id="Barni+ (2024)">[Barni+ (2024)]</a>    
  G. Barni, S. Blasi and M. Vanvlasselaer,  
  *The hydrodynamics of inverse phase transitions*,  
  [JCAP **10** (2024), 042](https://doi.org/10.1088/1475-7516/2024/10/042) 
  [[arXiv:2406.01596 [hep-ph]](https://arxiv.org/abs/2406.01596)].

<a id="Barni+ (2025a)">[Barni+ (2025a)]</a>  
  G. Barni, S. Blasi and M. Vanvlasselaer,  
  *Inverse bubbles from broken symmetries*,  
  [[arXiv:2503.01951 [hep-ph]](https://arxiv.org/abs/2503.01951)].

<a id="Barni+ (2025b)">[Barni+ (2025b)]</a>  
  G. Barni, S. Blasi, E. Madge and M. Vanvlasselaer,  
  *Gravitational waves from the sound shell model: direct and inverse phase
  transitions in the early Universe*,  
  [[arXiv:2510.21439 [hep-ph]](https://arxiv.org/abs/2510.21439)].