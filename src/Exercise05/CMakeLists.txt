add_executable(05_1 1.cpp structs.hpp)
target_link_libraries(05_1 PUBLIC CONAN_PKG::rapidcsv CONAN_PKG::cxxopts project_config ariel_random lsn_libs project_warnings)

install(TARGETS 05_1 RUNTIME DESTINATION bin)
