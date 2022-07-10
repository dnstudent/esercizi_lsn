add_executable(01_1 1.cpp)
target_link_libraries(01_1 PUBLIC CONAN_PKG::rapidcsv CONAN_PKG::cxxopts project_config ariel_random lsn_libs project_warnings)

add_executable(01_2 2.cpp)
target_link_libraries(01_2 PUBLIC CONAN_PKG::rapidcsv CONAN_PKG::cxxopts project_config ariel_random lsn_libs project_warnings)

add_executable(01_3 3.cpp)
target_link_libraries(01_3 PUBLIC CONAN_PKG::rapidcsv CONAN_PKG::cxxopts project_config ariel_random lsn_libs project_warnings)

install(TARGETS 01_1 01_2 01_3 RUNTIME DESTINATION bin)

