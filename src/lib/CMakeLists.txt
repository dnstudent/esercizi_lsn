add_subdirectory(ariel_random)
add_subdirectory(genetic)

add_library(lsn_libs INTERFACE)
target_include_directories(lsn_libs INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(lsn_libs INTERFACE CONAN_PKG::rapidcsv)