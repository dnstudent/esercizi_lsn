add_executable(02_1 1.cpp)
target_link_libraries(02_1 PUBLIC CONAN_PKG::rapidcsv CONAN_PKG::cxxopts project_config ariel_random lsn_libs project_warnings)

add_executable(02_2 2.cpp)
target_link_libraries(02_2 PUBLIC CONAN_PKG::rapidcsv CONAN_PKG::cxxopts CONAN_PKG::indicators project_config ariel_random lsn_libs project_warnings)

install(TARGETS 02_1 02_2 RUNTIME DESTINATION bin)