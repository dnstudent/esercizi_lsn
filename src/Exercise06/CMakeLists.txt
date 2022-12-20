add_executable(06_warmup equilibrate.cpp)
target_link_libraries(06_warmup PRIVATE CONAN_PKG::rapidcsv CONAN_PKG::cxxopts CONAN_PKG::indicators project_config ariel_random lsn_libs project_warnings)

add_executable(06_correlation correlation.cpp)
target_link_libraries(06_correlation PRIVATE CONAN_PKG::rapidcsv CONAN_PKG::cxxopts project_config lsn_libs project_warnings)

add_executable(06_measure measure.cpp)
target_link_libraries(06_measure PRIVATE CONAN_PKG::rapidcsv CONAN_PKG::cxxopts CONAN_PKG::indicators project_config ariel_random lsn_libs project_warnings)

install(TARGETS 06_warmup 06_correlation 06_measure RUNTIME DESTINATION bin)