add_subdirectory(ariel_random)

add_library(lsn_libs INTERFACE estimators/variance.hpp)
target_include_directories(lsn_libs INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})