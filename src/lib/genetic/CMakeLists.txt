add_library(genetic INTERFACE)
target_link_libraries(genetic INTERFACE project_warnings CONAN_PKG::indicators ${MPI_TARGETS})
target_include_directories(genetic INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})