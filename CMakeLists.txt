cmake_minimum_required(VERSION 3.0.2)

project(boundstates)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# ADD_DEFINITIONS(-DBOOST_ALL_DYN_LINK)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -lboost_system -g -lpthread -pthread")
set(COMPTOOLS_INCLUDE_DIR "~/software/comp-tools/")
set(COMPTOOLS_CPP
	"${COMPTOOLS_INCLUDE_DIR}comptools/function.cpp"
	"${COMPTOOLS_INCLUDE_DIR}comptools/math.cpp"
)

find_package(Boost 1.62 COMPONENTS system filesystem program_options REQUIRED)
find_package(LAPACK REQUIRED)

set(LCP_CPP
	"src/lcp.cpp"
	"src/lcp-app.cpp"
	"src/schrodinger-box.cpp"
)
add_executable(lcp ${LCP_CPP} ${COMPTOOLS_CPP})
target_include_directories(lcp PUBLIC
	include
	"${COMPTOOLS_INCLUDE_DIR}"
)
target_link_libraries(lcp
	"${LAPACK_LIBRARIES}"
)

set(DVR_CPP
	"src/dvr.cpp"
	"src/boundstates-app.cpp"
	"src/schrodinger-box.cpp"
)

add_executable(dvr ${DVR_CPP} ${COMPTOOLS_CPP})
target_include_directories(dvr PUBLIC
	include
	"${COMPTOOLS_INCLUDE_DIR}"
    "${Boost_INCLUDE_DIR}"
)
target_link_libraries(dvr
	"${Boost_LIBRARIES}"
	"${LAPACK_LIBRARIES}"
)

