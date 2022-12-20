add_executable(03_1 1.cpp)
target_link_libraries(03_1 PUBLIC CONAN_PKG::rapidcsv CONAN_PKG::cxxopts project_config ariel_random lsn_libs project_warnings)

install(TARGETS 03_1 RUNTIME DESTINATION bin)