cmake_minimum_required(VERSION 3.0.2)

project(boundstates)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

ADD_DEFINITIONS(-DBOOST_ALL_DYN_LINK)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -lboost_system -g")
set(COMPTOOLS_INCLUDE_DIR "~/software/comp-tools/")

find_package(Boost 1.62 COMPONENTS system filesystem program_options REQUIRED)
find_package(LAPACK REQUIRED)

include_directories(
	include
	"${COMPTOOLS_INCLUDE_DIR}"
    "${Boost_INCLUDE_DIR}"
)

file(GLOB SOURCES "src/*.cpp")

add_executable(boundstates ${SOURCES} "${COMPTOOLS_INCLUDE_DIR}comptools/function.cpp" "${COMPTOOLS_INCLUDE_DIR}comptools/math.cpp")

target_link_libraries(boundstates
	"${Boost_LIBRARIES}"
	"${LAPACK_LIBRARIES}"
)
