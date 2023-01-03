########### AGGREGATED COMPONENTS AND DEPENDENCIES FOR THE MULTI CONFIG #####################
#############################################################################################

set(indicators_COMPONENT_NAMES "")
set(indicators_FIND_DEPENDENCY_NAMES "")

########### VARIABLES #######################################################################
#############################################################################################
set(indicators_PACKAGE_FOLDER_RELEASE "/Users/davidenicoli/.conan/data/indicators/2.2/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9")
set(indicators_BUILD_MODULES_PATHS_RELEASE )


set(indicators_INCLUDE_DIRS_RELEASE "${indicators_PACKAGE_FOLDER_RELEASE}/include")
set(indicators_RES_DIRS_RELEASE "${indicators_PACKAGE_FOLDER_RELEASE}/res")
set(indicators_DEFINITIONS_RELEASE )
set(indicators_SHARED_LINK_FLAGS_RELEASE )
set(indicators_EXE_LINK_FLAGS_RELEASE )
set(indicators_OBJECTS_RELEASE )
set(indicators_COMPILE_DEFINITIONS_RELEASE )
set(indicators_COMPILE_OPTIONS_C_RELEASE )
set(indicators_COMPILE_OPTIONS_CXX_RELEASE )
set(indicators_LIB_DIRS_RELEASE "${indicators_PACKAGE_FOLDER_RELEASE}/lib")
set(indicators_LIBS_RELEASE )
set(indicators_SYSTEM_LIBS_RELEASE )
set(indicators_FRAMEWORK_DIRS_RELEASE "${indicators_PACKAGE_FOLDER_RELEASE}/Frameworks")
set(indicators_FRAMEWORKS_RELEASE )
set(indicators_BUILD_DIRS_RELEASE "${indicators_PACKAGE_FOLDER_RELEASE}/")

# COMPOUND VARIABLES
set(indicators_COMPILE_OPTIONS_RELEASE
    "$<$<COMPILE_LANGUAGE:CXX>:${indicators_COMPILE_OPTIONS_CXX_RELEASE}>"
    "$<$<COMPILE_LANGUAGE:C>:${indicators_COMPILE_OPTIONS_C_RELEASE}>")
set(indicators_LINKER_FLAGS_RELEASE
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,SHARED_LIBRARY>:${indicators_SHARED_LINK_FLAGS_RELEASE}>"
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,MODULE_LIBRARY>:${indicators_SHARED_LINK_FLAGS_RELEASE}>"
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,EXECUTABLE>:${indicators_EXE_LINK_FLAGS_RELEASE}>")


set(indicators_COMPONENTS_RELEASE )