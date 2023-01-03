

# Conan automatically generated toolchain file
# DO NOT EDIT MANUALLY, it will be overwritten

# Avoid including toolchain file several times (bad if appending to variables like
#   CMAKE_CXX_FLAGS. See https://github.com/android/ndk/issues/323
include_guard()

message(STATUS "Using Conan toolchain: ${CMAKE_CURRENT_LIST_FILE}")

if(${CMAKE_VERSION} VERSION_LESS "3.15")
    message(FATAL_ERROR "The 'CMakeToolchain' generator only works with CMake >= 3.15")
endif()










# Set the architectures for which to build.
set(CMAKE_OSX_ARCHITECTURES arm64 CACHE STRING "" FORCE)
# Setting CMAKE_OSX_SYSROOT SDK, when using Xcode generator the name is enough
# but full path is necessary for others
set(CMAKE_OSX_SYSROOT macosx CACHE STRING "" FORCE)
set(BITCODE "")
set(FOBJC_ARC "")
set(VISIBILITY "")
#Check if Xcode generator is used, since that will handle these flags automagically
if(CMAKE_GENERATOR MATCHES "Xcode")
  message(DEBUG "Not setting any manual command-line buildflags, since Xcode is selected as generator.")
else()
    string(APPEND CONAN_C_FLAGS " ${BITCODE} ${FOBJC_ARC}")
    string(APPEND CONAN_CXX_FLAGS " ${BITCODE} ${VISIBILITY} ${FOBJC_ARC}")
endif()

string(APPEND CONAN_CXX_FLAGS " -stdlib=libc++")


message(STATUS "Conan toolchain: C++ Standard 17 with extensions OFF")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Extra c, cxx, linkflags and defines


if(DEFINED CONAN_CXX_FLAGS)
  string(APPEND CMAKE_CXX_FLAGS_INIT " ${CONAN_CXX_FLAGS}")
endif()
if(DEFINED CONAN_C_FLAGS)
  string(APPEND CMAKE_C_FLAGS_INIT " ${CONAN_C_FLAGS}")
endif()
if(DEFINED CONAN_SHARED_LINKER_FLAGS)
  string(APPEND CMAKE_SHARED_LINKER_FLAGS_INIT " ${CONAN_SHARED_LINKER_FLAGS}")
endif()
if(DEFINED CONAN_EXE_LINKER_FLAGS)
  string(APPEND CMAKE_EXE_LINKER_FLAGS_INIT " ${CONAN_EXE_LINKER_FLAGS}")
endif()

get_property( _CMAKE_IN_TRY_COMPILE GLOBAL PROPERTY IN_TRY_COMPILE )
if(_CMAKE_IN_TRY_COMPILE)
    message(STATUS "Running toolchain IN_TRY_COMPILE")
    return()
endif()

set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ON)

# Definition of CMAKE_MODULE_PATH
# Explicitly defined "builddirs" of "host" dependencies
list(PREPEND CMAKE_MODULE_PATH "/Users/davidenicoli/.conan/data/catch2/3.1.0/_/_/package/63312d5711219c8db2043ea7a697bcde8b3d4e92/lib/cmake/Catch2")
# The root (which is the default builddirs) path of dependencies in the host context
list(PREPEND CMAKE_MODULE_PATH "/Users/davidenicoli/.conan/data/cxxopts/3.0.0/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9/" "/Users/davidenicoli/.conan/data/rapidcsv/8.64/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9/" "/Users/davidenicoli/.conan/data/indicators/2.2/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9/")
# the generators folder (where conan generates files, like this toolchain)
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

# Definition of CMAKE_PREFIX_PATH, CMAKE_XXXXX_PATH
# The explicitly defined "builddirs" of "host" context dependencies must be in PREFIX_PATH
list(PREPEND CMAKE_PREFIX_PATH "/Users/davidenicoli/.conan/data/catch2/3.1.0/_/_/package/63312d5711219c8db2043ea7a697bcde8b3d4e92/lib/cmake/Catch2")
# The Conan local "generators" folder, where this toolchain is saved.
list(PREPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_LIST_DIR} )
list(PREPEND CMAKE_LIBRARY_PATH "/Users/davidenicoli/.conan/data/cxxopts/3.0.0/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9/lib" "/Users/davidenicoli/.conan/data/rapidcsv/8.64/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9/lib" "/Users/davidenicoli/.conan/data/catch2/3.1.0/_/_/package/63312d5711219c8db2043ea7a697bcde8b3d4e92/lib" "/Users/davidenicoli/.conan/data/indicators/2.2/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9/lib")
list(PREPEND CMAKE_FRAMEWORK_PATH "/Users/davidenicoli/.conan/data/cxxopts/3.0.0/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9/Frameworks" "/Users/davidenicoli/.conan/data/rapidcsv/8.64/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9/Frameworks" "/Users/davidenicoli/.conan/data/indicators/2.2/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9/Frameworks")
list(PREPEND CMAKE_INCLUDE_PATH "/Users/davidenicoli/.conan/data/cxxopts/3.0.0/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9/include" "/Users/davidenicoli/.conan/data/rapidcsv/8.64/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9/include" "/Users/davidenicoli/.conan/data/catch2/3.1.0/_/_/package/63312d5711219c8db2043ea7a697bcde8b3d4e92/include" "/Users/davidenicoli/.conan/data/indicators/2.2/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9/include")



if (DEFINED ENV{PKG_CONFIG_PATH})
set(ENV{PKG_CONFIG_PATH} "/Users/davidenicoli/Local_Workspace/Uni/LabSiNum/esercizi_lsn/build/generators:$ENV{PKG_CONFIG_PATH}")
else()
set(ENV{PKG_CONFIG_PATH} "/Users/davidenicoli/Local_Workspace/Uni/LabSiNum/esercizi_lsn/build/generators:")
endif()




# Variables
# Variables  per configuration


# Preprocessor definitions
# Preprocessor definitions per configuration
