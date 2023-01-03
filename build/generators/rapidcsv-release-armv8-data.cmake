########### AGGREGATED COMPONENTS AND DEPENDENCIES FOR THE MULTI CONFIG #####################
#############################################################################################

set(rapidcsv_COMPONENT_NAMES "")
set(rapidcsv_FIND_DEPENDENCY_NAMES "")

########### VARIABLES #######################################################################
#############################################################################################
set(rapidcsv_PACKAGE_FOLDER_RELEASE "/Users/davidenicoli/.conan/data/rapidcsv/8.64/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9")
set(rapidcsv_BUILD_MODULES_PATHS_RELEASE )


set(rapidcsv_INCLUDE_DIRS_RELEASE "${rapidcsv_PACKAGE_FOLDER_RELEASE}/include")
set(rapidcsv_RES_DIRS_RELEASE "${rapidcsv_PACKAGE_FOLDER_RELEASE}/res")
set(rapidcsv_DEFINITIONS_RELEASE )
set(rapidcsv_SHARED_LINK_FLAGS_RELEASE )
set(rapidcsv_EXE_LINK_FLAGS_RELEASE )
set(rapidcsv_OBJECTS_RELEASE )
set(rapidcsv_COMPILE_DEFINITIONS_RELEASE )
set(rapidcsv_COMPILE_OPTIONS_C_RELEASE )
set(rapidcsv_COMPILE_OPTIONS_CXX_RELEASE )
set(rapidcsv_LIB_DIRS_RELEASE "${rapidcsv_PACKAGE_FOLDER_RELEASE}/lib")
set(rapidcsv_LIBS_RELEASE )
set(rapidcsv_SYSTEM_LIBS_RELEASE )
set(rapidcsv_FRAMEWORK_DIRS_RELEASE "${rapidcsv_PACKAGE_FOLDER_RELEASE}/Frameworks")
set(rapidcsv_FRAMEWORKS_RELEASE )
set(rapidcsv_BUILD_DIRS_RELEASE "${rapidcsv_PACKAGE_FOLDER_RELEASE}/")

# COMPOUND VARIABLES
set(rapidcsv_COMPILE_OPTIONS_RELEASE
    "$<$<COMPILE_LANGUAGE:CXX>:${rapidcsv_COMPILE_OPTIONS_CXX_RELEASE}>"
    "$<$<COMPILE_LANGUAGE:C>:${rapidcsv_COMPILE_OPTIONS_C_RELEASE}>")
set(rapidcsv_LINKER_FLAGS_RELEASE
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,SHARED_LIBRARY>:${rapidcsv_SHARED_LINK_FLAGS_RELEASE}>"
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,MODULE_LIBRARY>:${rapidcsv_SHARED_LINK_FLAGS_RELEASE}>"
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,EXECUTABLE>:${rapidcsv_EXE_LINK_FLAGS_RELEASE}>")


set(rapidcsv_COMPONENTS_RELEASE )