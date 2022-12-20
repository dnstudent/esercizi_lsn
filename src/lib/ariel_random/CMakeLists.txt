add_library(ariel_random STATIC ariel_random/ariel_random.cpp)
target_include_directories(ariel_random
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(ariel_random PROPERTIES CXX_EXTENSIONS OFF)
install(TARGETS ariel_random LIBRARY DESTINATION lib)
