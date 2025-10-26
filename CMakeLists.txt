cmake_minimum_required(VERSION 3.24...3.31)

project(
    inverse_pt_calculator 
    VERSION 
        0.1
    DESCRIPTION 
        "Compute fluid profiles and the gravitational wave spectrum for (direct and) inverse phase transitions."
    LANGUAGES 
        CXX
)

# use Release as default build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Build config
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MACOSX_RPATH ON)

# inverse_pt_calculator build options
option(BUILD_PYTHONLIB "Build the Python interface" ON)
option(BUILD_DOCS "Build Doxygen documentation" ON)
option(BUILD_TEST "Build the test executable" ON)

# cmake compiler options
option(CMAKE_POSITION_INDEPENDENT_CODE "Enable PIC" ON) # -fPIC
option(CMAKE_INTERPROCEDURAL_OPTIMIZATION "Enable IPO/LTO" ON) # -flto
# compiler options
option(OPTIMIZE_NATIVE "Optimize code for the machine architecture (-march=native)" ON)
option(DEVELOPMENT_WARNINGS "Switch on warnings for development (-Wall)" OFF)
option(VECTORIZATION_REPORTS "Enable loop vectorization reports (-Rpass)" OFF)
# option for reducing doxygen output
option(QUIET_DOXYGEN "Only display warnings and errors when creating the documentation" ON)

# Required dependencies
find_package(Boost CONFIG REQUIRED COMPONENTS log)
find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)
# Optional dependencies
find_package(Python3 COMPONENTS Interpreter Development.Module)
find_package(pybind11)
if(NOT (Python3_FOUND AND pybind11_FOUND))
    if(BUILD_PYTHONLIB)
        message(FATAL_ERROR
            "Python3 and/or pybind11 not found."
            " Please either provide help for locating the packages (e.g. Python3_ROOT_DIR and pybind11_DIR)"
            " or set BUILD_PYTHONLIB=OFF."
        )
    else()
        message(STATUS
            "Python3 and/or pybind11 not found."
            " You will not be able to build the Python interface."
        )
    endif()
endif()
find_package(Doxygen)
if(NOT Doxygen_FOUND)
    if(BUILD_DOCS)
        message(FATAL_ERROR
            "Doxygen not found."
            " Please either provide help for locating the package or set BUILD_DOCS=OFF."
        )
    else()
        message(STATUS
            "Doxygen not found. You will not be able to build the documentation."
        )
    endif()
endif()

# set compiler options
add_library(inverse_pt_compiler_flags INTERFACE)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    if (OPTIMIZE_NATIVE)
        message(STATUS "Code is optimized for the machine-specific architecture.")
        target_compile_options(inverse_pt_compiler_flags INTERFACE -march=native)
    endif()
    if (DEVELOPMENT_WARNINGS)
        message(STATUS "Compilation warnings are enabled.")
        target_compile_options(inverse_pt_compiler_flags INTERFACE -Wall)      
    endif()
    if (VECTORIZATION_REPORTS)
        message(STATUS "Vectorization reports are enabled.")
        target_compile_options(inverse_pt_compiler_flags INTERFACE
            -Rpass=loop-vectorize
            -Rpass-missed=loop-vectorize
        )
    endif()
else()
    message(
        WARNING 
        "Compiler '${CMAKE_CXX_COMPILER_ID}' is not explicitly supported (supported compilers are GNU and Clang).\n"
        "The compilation flag options DEVELOPMENT_WARNINGS, OPTIMIZE_NATIVE, and VECTORIZATION_REPORTS are ignored.\n"
        "If you want to enable similar flags, please pass them manually using:\n"
        "    -DCMAKE_CXX_FLAGS=\"<your flags>\"\n"
        "or set environment variables like:\n"
        "    CXXFLAGS=\"<your flags>\"\n"
    )
endif()

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test)
# Python site-packages directory (relative to lib)
set(PYSITE_SUBDIRECTORY python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages)

# define variable for setting relative RPATH
if(APPLE)
    set(RPATH_EXECUTABLE "@executable_path")
    set(RPATH_LIBRARY "@loader_path")
elseif(UNIX)
    set(RPATH_EXECUTABLE "\\$ORIGIN")
    set(RPATH_LIBRARY "\\$ORIGIN")
endif()

# C++ shared library
add_library(inverse_pt_calculator SHARED
    src/fluid_profile.cpp
    src/profile_calculator.cpp
    src/sound_shell_model.cpp
    src/sound_shell_spectrum.cpp
    src/utils.cpp
    src/settings.cpp
)
target_include_directories(inverse_pt_calculator 
    PUBLIC ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(inverse_pt_calculator
    PUBLIC Boost::log
    PUBLIC GSL::gsl GSL::gslcblas
    PUBLIC OpenMP::OpenMP_CXX
    PRIVATE inverse_pt_compiler_flags
)
# Install shared library + headers
install(TARGETS inverse_pt_calculator
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)

# Python interface
if(Python3_FOUND AND pybind11_FOUND)
    pybind11_add_module(pythonlib
        src/pythonlib/fluid_profile_python.cpp
        src/pythonlib/profile_calculator_python.cpp
        src/pythonlib/sound_shell_model_python.cpp
        src/pythonlib/sound_shell_spectrum_python.cpp
        src/pythonlib/utils_python.cpp
        src/pythonlib/settings_python.cpp
        src/pythonlib/inverse_pt_calculator.cpp
    )
    target_include_directories(pythonlib 
        PRIVATE ${PROJECT_SOURCE_DIR}/src/pythonlib
    )
    target_link_libraries(pythonlib 
        PRIVATE inverse_pt_calculator inverse_pt_compiler_flags
    )
    set_target_properties(pythonlib PROPERTIES
        OUTPUT_NAME "inverse_pt_calculator"
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${PYSITE_SUBDIRECTORY}
        BUILD_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY};${RPATH_EXECUTABLE}/../.."
        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${RPATH_EXECUTABLE}/../.."
    )
    message(STATUS "BUILD_PYTHONLIB for building the Python interface is set to ${BUILD_PYTHONLIB}.")
    if(BUILD_PYTHONLIB)
        # install
        message(STATUS "The Python interface will be built and installed to ${CMAKE_INSTALL_PREFIX}/lib/${PYSITE_SUBDIRECTORY}.")
        install(TARGETS pythonlib
            LIBRARY DESTINATION lib/${PYSITE_SUBDIRECTORY}
        )
    else()
        # do not build pythonlib by default unless requested
        message(STATUS "The Python inteface will not be built by default and will not be installed.")
        set_property(TARGET pythonlib PROPERTY EXCLUDE_FROM_ALL TRUE)
    endif()
else()
    message(STATUS "Python or pybind were not found. The Python interface cannot be built.")
endif()

# Test executable
add_executable(test test/test.cpp)
target_link_libraries(test
    PRIVATE inverse_pt_calculator inverse_pt_compiler_flags
)
set_target_properties(test PROPERTIES
    BUILD_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY};${RPATH_EXECUTABLE}/../lib"
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${RPATH_EXECUTABLE}/../lib"
)
message(STATUS "BUILD_TEST for building the test executable is set to ${BUILD_TEST}")
if(BUILD_TEST)
    message(STATUS "The test executable will be built and installed to ${CMAKE_INSTALL_PREFIX}/bin.")
    install(TARGETS test DESTINATION bin)
else()
    # do not build test by default unless requested
    message(STATUS "The test executable will not be built by default and will not be installed.")
    set_property(TARGET test PROPERTY EXCLUDE_FROM_ALL TRUE)
endif()

# Documentation
if(DOXYGEN_FOUND)
    # prepare config file
    set(DOXYGEN_CONFIG_OUT ${CMAKE_BINARY_DIR}/doxygen.config)
    file(READ ${PROJECT_SOURCE_DIR}/doxygen.config DOXYGEN_CONFIG_CONTENTS)
    file(WRITE ${DOXYGEN_CONFIG_OUT} "${DOXYGEN_CONFIG_CONTENTS}\n")
    file(APPEND ${DOXYGEN_CONFIG_OUT} "OUTPUT_DIRECTORY = ${CMAKE_BINARY_DIR}/docs\n")
    file(APPEND ${DOXYGEN_CONFIG_OUT} "INPUT = ${PROJECT_SOURCE_DIR}/include/inverse_pt ${PROJECT_SOURCE_DIR}/src\n")
    if(QUIET_DOXYGEN)
        file(APPEND ${DOXYGEN_CONFIG_OUT} "QUIET = YES\n")
    endif()
    # Doxygen docs target
    add_custom_target(docs
        COMMAND  Doxygen::doxygen ${DOXYGEN_CONFIG_OUT}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Generating API documentation into build/docs"
    )
    message(STATUS "BUILD_DOCS for building the documentation is set to ${BUILD_DOCS}")
    # build docs by default (and install them) when requested
    if(BUILD_DOCS)
        message(STATUS "The documentation will be built and installed to ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/.")
        set_property(TARGET docs PROPERTY EXCLUDE_FROM_ALL FALSE)
        install(DIRECTORY ${CMAKE_BINARY_DIR}/docs DESTINATION share/${PROJECT_NAME}/)
    else()
        message(STATUS "The documentation will not be built by default and will not be installed.")
    endif()
endif()

# target to build everything
add_custom_target(everything DEPENDS inverse_pt_calculator test)
if(Python3_FOUND AND pybind11_FOUND)
    add_dependencies(everything pythonlib)
endif()
if(DOXYGEN_FOUND)
    add_dependencies(everything docs)
endif()
