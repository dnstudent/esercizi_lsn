# Avoid multiple calls to find_package to append duplicated properties to the targets
include_guard()########### VARIABLES #######################################################################
#############################################################################################
set(rapidcsv_FRAMEWORKS_FOUND_RELEASE "") # Will be filled later
conan_find_apple_frameworks(rapidcsv_FRAMEWORKS_FOUND_RELEASE "${rapidcsv_FRAMEWORKS_RELEASE}" "${rapidcsv_FRAMEWORK_DIRS_RELEASE}")

set(rapidcsv_LIBRARIES_TARGETS "") # Will be filled later


######## Create an interface target to contain all the dependencies (frameworks, system and conan deps)
if(NOT TARGET rapidcsv_DEPS_TARGET)
    add_library(rapidcsv_DEPS_TARGET INTERFACE IMPORTED)
endif()

set_property(TARGET rapidcsv_DEPS_TARGET
             PROPERTY INTERFACE_LINK_LIBRARIES
             $<$<CONFIG:Release>:${rapidcsv_FRAMEWORKS_FOUND_RELEASE}>
             $<$<CONFIG:Release>:${rapidcsv_SYSTEM_LIBS_RELEASE}>
             $<$<CONFIG:Release>:>
             APPEND)

####### Find the libraries declared in cpp_info.libs, create an IMPORTED target for each one and link the
####### rapidcsv_DEPS_TARGET to all of them
conan_package_library_targets("${rapidcsv_LIBS_RELEASE}"    # libraries
                              "${rapidcsv_LIB_DIRS_RELEASE}" # package_libdir
                              rapidcsv_DEPS_TARGET
                              rapidcsv_LIBRARIES_TARGETS  # out_libraries_targets
                              "_RELEASE"
                              "rapidcsv")    # package_name

# FIXME: What is the result of this for multi-config? All configs adding themselves to path?
set(CMAKE_MODULE_PATH ${rapidcsv_BUILD_DIRS_RELEASE} ${CMAKE_MODULE_PATH})


########## GLOBAL TARGET PROPERTIES Release ########################################
    set_property(TARGET rapidcsv::rapidcsv
                 PROPERTY INTERFACE_LINK_LIBRARIES
                 $<$<CONFIG:Release>:${rapidcsv_OBJECTS_RELEASE}>
                 $<$<CONFIG:Release>:${rapidcsv_LIBRARIES_TARGETS}>
                 APPEND)

    if("${rapidcsv_LIBS_RELEASE}" STREQUAL "")
        # If the package is not declaring any "cpp_info.libs" the package deps, system libs,
        # frameworks etc are not linked to the imported targets and we need to do it to the
        # global target
        set_property(TARGET rapidcsv::rapidcsv
                     PROPERTY INTERFACE_LINK_LIBRARIES
                     rapidcsv_DEPS_TARGET
                     APPEND)
    endif()

    set_property(TARGET rapidcsv::rapidcsv
                 PROPERTY INTERFACE_LINK_OPTIONS
                 $<$<CONFIG:Release>:${rapidcsv_LINKER_FLAGS_RELEASE}> APPEND)
    set_property(TARGET rapidcsv::rapidcsv
                 PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                 $<$<CONFIG:Release>:${rapidcsv_INCLUDE_DIRS_RELEASE}> APPEND)
    set_property(TARGET rapidcsv::rapidcsv
                 PROPERTY INTERFACE_COMPILE_DEFINITIONS
                 $<$<CONFIG:Release>:${rapidcsv_COMPILE_DEFINITIONS_RELEASE}> APPEND)
    set_property(TARGET rapidcsv::rapidcsv
                 PROPERTY INTERFACE_COMPILE_OPTIONS
                 $<$<CONFIG:Release>:${rapidcsv_COMPILE_OPTIONS_RELEASE}> APPEND)

########## For the modules (FindXXX)
set(rapidcsv_LIBRARIES_RELEASE rapidcsv::rapidcsv)
