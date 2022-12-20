add_executable(10_2 2.cpp)
target_link_libraries(10_2 PRIVATE CONAN_PKG::cxxopts CONAN_PKG::rapidcsv CONAN_PKG::indicators project_config ariel_random lsn_libs genetic project_warnings)
target_compile_definitions(10_2 PRIVATE USE_MPI=1)

install(TARGETS 10_2 RUNTIME DESTINATION bin)