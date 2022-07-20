add_library(ariel_mdmc SHARED ariel_mdmc/MD_MC.cpp)

target_link_libraries(ariel_mdmc PRIVATE ariel_random)

target_include_directories(ariel_mdmc
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(ariel_mdmc PROPERTIES CXX_EXTENSIONS OFF)
install(TARGETS ariel_mdmc LIBRARY DESTINATION lib)
