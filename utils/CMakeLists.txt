# TODO remove once we increase the cmake requirement
if (CMAKE_VERSION VERSION_LESS 2.8.9)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -DPIC")
else()
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

set(atools_SRCS atools.cc)
add_executable(atools ${atools_SRCS})
target_link_libraries(atools utils ${Boost_LIBRARIES} ${ZLIB_LIBRARIES})

set(utils_STAT_SRCS
    alignment_io.h
    array2d.h
    basic_exception_type.h
    exception.h
    exception_type.h
    filelib.h
    gzstream.h
    null_deleter.h
    optional.h
    alignment_io.cc
    basic_exception_type.cc
    exception_type.cc
    filelib.cc
    gzstream.cc)

add_library(utils STATIC ${utils_STAT_SRCS})
