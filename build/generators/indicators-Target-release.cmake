# Avoid multiple calls to find_package to append duplicated properties to the targets
include_guard()########### VARIABLES #######################################################################
#############################################################################################
set(indicators_FRAMEWORKS_FOUND_RELEASE "") # Will be filled later
conan_find_apple_frameworks(indicators_FRAMEWORKS_FOUND_RELEASE "${indicators_FRAMEWORKS_RELEASE}" "${indicators_FRAMEWORK_DIRS_RELEASE}")

set(indicators_LIBRARIES_TARGETS "") # Will be filled later


######## Create an interface target to contain all the dependencies (frameworks, system and conan deps)
if(NOT TARGET indicators_DEPS_TARGET)
    add_library(indicators_DEPS_TARGET INTERFACE IMPORTED)
endif()

set_property(TARGET indicators_DEPS_TARGET
             PROPERTY INTERFACE_LINK_LIBRARIES
             $<$<CONFIG:Release>:${indicators_FRAMEWORKS_FOUND_RELEASE}>
             $<$<CONFIG:Release>:${indicators_SYSTEM_LIBS_RELEASE}>
             $<$<CONFIG:Release>:>
             APPEND)

####### Find the libraries declared in cpp_info.libs, create an IMPORTED target for each one and link the
####### indicators_DEPS_TARGET to all of them
conan_package_library_targets("${indicators_LIBS_RELEASE}"    # libraries
                              "${indicators_LIB_DIRS_RELEASE}" # package_libdir
                              indicators_DEPS_TARGET
                              indicators_LIBRARIES_TARGETS  # out_libraries_targets
                              "_RELEASE"
                              "indicators")    # package_name

# FIXME: What is the result of this for multi-config? All configs adding themselves to path?
set(CMAKE_MODULE_PATH ${indicators_BUILD_DIRS_RELEASE} ${CMAKE_MODULE_PATH})


########## GLOBAL TARGET PROPERTIES Release ########################################
    set_property(TARGET indicators::indicators
                 PROPERTY INTERFACE_LINK_LIBRARIES
                 $<$<CONFIG:Release>:${indicators_OBJECTS_RELEASE}>
                 $<$<CONFIG:Release>:${indicators_LIBRARIES_TARGETS}>
                 APPEND)

    if("${indicators_LIBS_RELEASE}" STREQUAL "")
        # If the package is not declaring any "cpp_info.libs" the package deps, system libs,
        # frameworks etc are not linked to the imported targets and we need to do it to the
        # global target
        set_property(TARGET indicators::indicators
                     PROPERTY INTERFACE_LINK_LIBRARIES
                     indicators_DEPS_TARGET
                     APPEND)
    endif()

    set_property(TARGET indicators::indicators
                 PROPERTY INTERFACE_LINK_OPTIONS
                 $<$<CONFIG:Release>:${indicators_LINKER_FLAGS_RELEASE}> APPEND)
    set_property(TARGET indicators::indicators
                 PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                 $<$<CONFIG:Release>:${indicators_INCLUDE_DIRS_RELEASE}> APPEND)
    set_property(TARGET indicators::indicators
                 PROPERTY INTERFACE_COMPILE_DEFINITIONS
                 $<$<CONFIG:Release>:${indicators_COMPILE_DEFINITIONS_RELEASE}> APPEND)
    set_property(TARGET indicators::indicators
                 PROPERTY INTERFACE_COMPILE_OPTIONS
                 $<$<CONFIG:Release>:${indicators_COMPILE_OPTIONS_RELEASE}> APPEND)

########## For the modules (FindXXX)
set(indicators_LIBRARIES_RELEASE indicators::indicators)
