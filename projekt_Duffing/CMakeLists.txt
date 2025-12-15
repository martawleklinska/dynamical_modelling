cmake_minimum_required(VERSION 3.10)
project(ODE_Templates)

set(CMAKE_CXX_STANDARD 17)

include_directories(${PROJECT_SOURCE_DIR}/include)

find_package(GSL REQUIRED)

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(FFTW3 IMPORTED_TARGET fftw3)
endif()

if(NOT FFTW3_FOUND)
    find_path(FFTW3_INCLUDE_DIR fftw3.h)
    find_library(FFTW3_LIBRARY NAMES fftw3)
    
    if(FFTW3_INCLUDE_DIR AND FFTW3_LIBRARY)
        set(FFTW3_FOUND TRUE)
        set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
        set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})
    endif()
endif()

if(NOT FFTW3_FOUND)
    message(FATAL_ERROR "FFTW3 not found. Please install FFTW3 development package.")
endif()

file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_executable(main ${SOURCES})

target_include_directories(main PRIVATE ${FFTW3_INCLUDE_DIRS})

target_link_libraries(main
    PRIVATE
        GSL::gsl
        GSL::gslcblas
        ${FFTW3_LIBRARIES}
)
