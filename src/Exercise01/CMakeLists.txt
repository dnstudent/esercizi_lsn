add_executable(01_1 1.1.cpp)
target_link_libraries(01_1 PUBLIC ${COMMON_LIBS})
target_include_directories(01_1 PUBLIC ${COMMON_INCLUDE})

add_executable(01_1_3 1.3.cpp)
target_link_libraries(01_1_3 PUBLIC ${COMMON_LIBS})
target_include_directories(01_1_3 PUBLIC ${COMMON_INCLUDE})

add_executable(01_2 2.cpp)
target_link_libraries(01_2 PUBLIC ${COMMON_LIBS})
target_include_directories(01_2 PUBLIC ${COMMON_INCLUDE})
