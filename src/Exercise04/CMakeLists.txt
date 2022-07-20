add_executable(04_1 1.cpp)
target_link_libraries(04_1 PUBLIC CONAN_PKG::rapidcsv CONAN_PKG::cxxopts project_config ariel_random lsn_libs project_warnings)

#install(TARGETS 04_1 RUNTIME DESTINATION bin)