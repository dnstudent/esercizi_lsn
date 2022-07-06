add_library(project_config INTERFACE)
target_include_directories(project_config INTERFACE "${PROJECT_SOURCE_DIR}/src/include")

#set(BUILD::warnings "project_options;project_warnings")
#add_library(project_devel INTERFACE project_options project_warnings)
#add_library(project::devel ALIAS project_devel)

configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/include/config.hpp.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/config.hpp")

add_subdirectory(extern)
add_subdirectory(lib)
add_subdirectory(Exercise01)