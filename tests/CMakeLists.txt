add_executable(tests TestEstimators.cpp)
target_link_libraries(tests PRIVATE CONAN_PKG::catch2 ariel_random lsn_libs)

include(CTest)
include(Catch)
catch_discover_tests(tests)