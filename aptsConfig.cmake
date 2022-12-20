# apts cmake module
#
# This module sets the target:
#
#   apts
#
# In addition, it sets the following variables:
#
#   apts_FOUND - true if apts found
#   apts_VERSION - apts's version
#   apts_INCLUDE_DIRS - the directory containing apts headers
#
# The following support targets are defined to simplify things:
#
#   apts::compiler_warnings - enable compiler warnings
#   apts::assert - enable apts assertions
#   apts::debug - enable all assertions (slow)

include(CMakeFindDependencyMacro)

# Define target "apts"

if(NOT TARGET apts)
    include("${CMAKE_CURRENT_LIST_DIR}/aptsTargets.cmake")
    get_target_property(apts_INCLUDE_DIRS apts INTERFACE_INCLUDE_DIRECTORIES)
endif()

# Find dependencies

find_dependency(xtensor)
find_dependency(prrng)

# Define support target "apts::compiler_warnings"

if(NOT TARGET apts::compiler_warnings)
    add_library(apts::compiler_warnings INTERFACE IMPORTED)
    if(MSVC)
        set_property(
            TARGET apts::compiler_warnings
            PROPERTY INTERFACE_COMPILE_OPTIONS
            /W4)
    else()
        set_property(
            TARGET apts::compiler_warnings
            PROPERTY INTERFACE_COMPILE_OPTIONS
            -Wall -Wextra -pedantic -Wno-unknown-pragmas)
    endif()
endif()

# Define support target "apts::assert"

if(NOT TARGET apts::assert)
    add_library(apts::assert INTERFACE IMPORTED)
    set_property(
        TARGET apts::assert
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        apts_ENABLE_ASSERT)
endif()

# Define support target "apts::debug"

if(NOT TARGET apts::debug)
    add_library(apts::debug INTERFACE IMPORTED)
    set_property(
        TARGET apts::debug
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        XTENSOR_ENABLE_ASSERT
        apts_ENABLE_ASSERT
        apts_ENABLE_DEBUG)
endif()
