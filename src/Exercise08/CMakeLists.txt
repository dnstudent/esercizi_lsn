add_executable(08_1 1.cpp)
target_link_libraries(08_1 PRIVATE CONAN_PKG::rapidcsv project_config ariel_random lsn_libs project_warnings)

add_executable(08_2 2.cpp)
target_link_libraries(08_2 PRIVATE CONAN_PKG::rapidcsv CONAN_PKG::cxxopts CONAN_PKG::indicators project_config ariel_random lsn_libs project_warnings)

add_executable(08_2_psi 2_psi.cpp)
target_link_libraries(08_2_psi PRIVATE CONAN_PKG::rapidcsv CONAN_PKG::cxxopts CONAN_PKG::indicators project_config ariel_random lsn_libs project_warnings)


install(TARGETS 08_1 08_2 08_2_psi RUNTIME DESTINATION bin)