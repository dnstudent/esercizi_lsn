# add_executable(04_1 1.cpp)
# target_link_libraries(04_1 PUBLIC CONAN_PKG::rapidcsv CONAN_PKG::cxxopts CONAN_PKG::indicators project_config ariel_random lsn_libs project_warnings)
add_executable(04_2 2.cpp)
target_link_libraries(04_2 PRIVATE CONAN_PKG::rapidcsv CONAN_PKG::cxxopts CONAN_PKG::indicators project_config ariel_random lsn_libs project_warnings)

install(TARGETS 04_2 RUNTIME DESTINATION bin)