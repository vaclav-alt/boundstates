cmake_minimum_required(VERSION 3.0.2)

project(boundstates)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -g -lboost_system -pthread")
set(COMPTOOLS_INCLUDE_DIR "~/software/comp-tools/include")

find_package(Boost 1.55 REQUIRED system log filesystem program_options)
find_package(Eigen3)

include_directories(
	include
	"${COMPTOOLS_INCLUDE_DIR}"
    "${Boost_INCLUDE_DIR}"
)

file(GLOB SOURCES "src/*.cpp")

add_executable(boundstates ${SOURCES})

target_link_libraries(boundstates
	"${Boost_LIBRARIES}"
    Eigen3::Eigen
)
