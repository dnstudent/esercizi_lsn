add_executable(tests GeneticTests.cpp EstimatorsTest.cpp AlgoTests.cpp VectorsTests.cpp MDTests.cpp RngTests.cpp TransitionsTests.cpp IsingTests.cpp MetaTests.cpp UtilsTests.cpp IntegratorTests.cpp Ex08Tests.cpp)
target_link_libraries(tests PRIVATE CONAN_PKG::catch2 genetic ariel_random lsn_libs project_config)

include(CTest)
include(Catch)
catch_discover_tests(tests)