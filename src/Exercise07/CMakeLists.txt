add_executable(07_2 2.cpp)
target_link_libraries(07_2 PRIVATE CONAN_PKG::rapidcsv CONAN_PKG::cxxopts project_config ariel_random lsn_libs project_warnings)

add_executable(07_4 4.cpp)
target_link_libraries(07_4 PRIVATE CONAN_PKG::rapidcsv CONAN_PKG::cxxopts CONAN_PKG::indicators project_config ariel_random lsn_libs project_warnings)


install(TARGETS 07_2 07_4 RUNTIME DESTINATION bin)