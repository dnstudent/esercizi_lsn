add_library(Ising SHARED ising.hpp)
set_target_properties(Ising PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(Ising PUBLIC Eigen3::Eigen)