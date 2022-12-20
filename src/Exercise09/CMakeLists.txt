add_executable(09_1 1.cpp)
target_link_libraries(09_1 PRIVATE CONAN_PKG::cxxopts CONAN_PKG::rapidcsv CONAN_PKG::indicators project_config ariel_random lsn_libs genetic project_warnings)
target_compile_definitions(09_1 PRIVATE USE_MPI=0)

install(TARGETS 09_1 RUNTIME DESTINATION bin)