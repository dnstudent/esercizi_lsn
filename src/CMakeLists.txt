add_library(project_config INTERFACE)
target_include_directories(project_config INTERFACE "${PROJECT_SOURCE_DIR}/src/include")

#set(BUILD::warnings "project_options;project_warnings")
#add_library(project_devel INTERFACE project_options project_warnings)
#add_library(project::devel ALIAS project_devel)

configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/include/config.hpp.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/config.hpp")

add_subdirectory(lib)
add_subdirectory(Exercise01)
add_subdirectory(Exercise02)
add_subdirectory(Exercise03)
add_subdirectory(Exercise04)
add_subdirectory(Exercise05)
add_subdirectory(Exercise06)
add_subdirectory(Exercise07)
add_subdirectory(Exercise08)
add_subdirectory(Exercise09)
add_subdirectory(Exercise10)