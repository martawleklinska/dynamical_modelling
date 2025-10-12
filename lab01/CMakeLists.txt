cmake_minimum_required(VERSION 3.10)
project(GSL_ODE_Solver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${PROJECT_SOURCE_DIR}/include)

find_package(GSL REQUIRED)

file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_executable(main ${SOURCES})

target_link_libraries(main PRIVATE GSL::gsl GSL::gslcblas)
