########## MACROS ###########################################################################
#############################################################################################

# Requires CMake > 3.15
if(${CMAKE_VERSION} VERSION_LESS "3.15")
    message(FATAL_ERROR "The 'CMakeDeps' generator only works with CMake >= 3.15")
endif()

include(${CMAKE_CURRENT_LIST_DIR}/cmakedeps_macros.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/rapidcsvTargets.cmake)
include(CMakeFindDependencyMacro)

check_build_type_defined()

foreach(_DEPENDENCY ${rapidcsv_FIND_DEPENDENCY_NAMES} )
    # Check that we have not already called a find_package with the transitive dependency
    if(NOT ${_DEPENDENCY}_FOUND)
        find_dependency(${_DEPENDENCY} REQUIRED ${${_DEPENDENCY}_FIND_MODE})
    endif()
endforeach()

set(rapidcsv_VERSION_STRING "8.64")
set(rapidcsv_INCLUDE_DIRS ${rapidcsv_INCLUDE_DIRS_RELEASE} )
set(rapidcsv_INCLUDE_DIR ${rapidcsv_INCLUDE_DIRS_RELEASE} )
set(rapidcsv_LIBRARIES ${rapidcsv_LIBRARIES_RELEASE} )
set(rapidcsv_DEFINITIONS ${rapidcsv_DEFINITIONS_RELEASE} )

# Only the first installed configuration is included to avoid the collision
foreach(_BUILD_MODULE ${rapidcsv_BUILD_MODULES_PATHS_RELEASE} )
    message(STATUS "Conan: Including build module from '${_BUILD_MODULE}'")
    include(${_BUILD_MODULE})
endforeach()


