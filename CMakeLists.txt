cmake_minimum_required(VERSION 2.8)
project(cdec)

add_definitions(-DKENLM_MAX_ORDER=6 -DHAVE_CONFIG_H)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_FLAGS "-Wall -std=c++11 -O3 -g")
set(METEOR_JAR "" CACHE FILEPATH "Path to meteor.jar")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# core packages
if(DEFINED ENV{BOOST_ROOT})
  set(Boost_NO_SYSTEM_PATHS ON)
endif()
set(Boost_REALPATH ON)

find_package(Boost COMPONENTS
                   regex
                   filesystem
                   serialization
                   program_options
                   unit_test_framework
                   system
                   thread REQUIRED)

include_directories(${Boost_INCLUDE_DIR})

# compression packages (primarily used by KenLM)
find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
  set(HAVE_ZLIB 1)
endif(ZLIB_FOUND)

# generate config.h
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)

add_subdirectory(utils)

set(CPACK_PACKAGE_VERSION_MAJOR "2015")
set(CPACK_PACKAGE_VERSION_MINOR "04")
set(CPACK_PACKAGE_VERSION_PATCH "26")
set(CPACK_SOURCE_GENERATOR "TBZ2")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_IGNORE_FILES
  "/.git/;/.gitignore;/Testing/;/build/;/.bzr/;~$;/CMakeCache.txt;/CMakeFiles/;${CPACK_SOURCE_IGNORE_FILES}")
include(CPack)
