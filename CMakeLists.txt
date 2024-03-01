cmake_minimum_required(VERSION 3.11)
project(class)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

# Options
option(ENABLE_TESTS "compile test executables" OFF)

# OpenMP
find_package(OpenMP)

# Sources
file( GLOB TOOLS_SRC
  ${PROJECT_SOURCE_DIR}/tools/*.c
)
file( GLOB SOURCE_SRC
  ${PROJECT_SOURCE_DIR}/source/*.c
)

set(CPP_SRC 
    cpp/ClassEngine.hh cpp/ClassEngine.cc 
    cpp/Engine.hh cpp/Engine.cc
)

set(TESTS
    test_loops
    test_loops_omp
    test_spectra
    test_transfer
    test_nonlinear
    test_perturbations
    test_thermodynamics
    test_background
    test_hyperspherical
)

file( GLOB EXTERNAL_SRC
    ${PROJECT_SOURCE_DIR}/external/RecfastCLASS/*.c
    ${PROJECT_SOURCE_DIR}/external/HyRec2020/*.c
    ${PROJECT_SOURCE_DIR}/external/heating/*.c
)
set(EXTERNAL_INCLUDE external/HyRec2020/)

## Main library
add_library(libclass ${TOOLS_SRC} ${SOURCE_SRC} ${RECFAST_SRC} ${EXTERNAL_SRC})
add_library(class::libclass ALIAS libclass)
set_target_properties(libclass PROPERTIES OUTPUT_NAME class)
set_target_properties(libclass PROPERTIES C_STANDARD 11)
set_property(TARGET libclass PROPERTY POSITION_INDEPENDENT_CODE ON)
# target_compile_options(libclass PRIVATE "-O4" "-ffast-math" "-g")
target_compile_definitions(libclass PUBLIC "__CLASSDIR__=\"${CMAKE_CURRENT_SOURCE_DIR}\"")
target_include_directories(libclass PUBLIC include ${EXTERNAL_INCLUDE} ${CMAKE_CURRENT_LIST_DIR}/external/RecfastCLASS ${CMAKE_CURRENT_LIST_DIR}/external/heating)
target_link_libraries(libclass PRIVATE m)
if(OpenMP_FOUND)
    target_link_libraries(libclass PRIVATE OpenMP::OpenMP_C)
endif()
if(ENABLE_HYREC)
    target_compile_definitions(libclass PUBLIC "HYREC")
    target_include_directories(libclass PUBLIC ${CMAKE_CURRENT_LIST_DIR}/external/HyRec2020/)
endif()

## Main Executable
add_executable(class main/class.c)
target_link_libraries(class PRIVATE libclass)

## Tests
if(ENABLE_TESTS)
    foreach(test ${TESTS})
        add_executable(${test} "test/${test}.c")
        # target_compile_options(${test} PRIVATE "-O4" "-ffast-math" "-g")
        target_link_libraries(${test} PRIVATE libclass)
        if(OpenMP_FOUND)
            target_link_libraries(${test} PRIVATE OpenMP::OpenMP_C)
        endif()
    endforeach()
endif()

## C++ library
add_library(libclass_cpp ${CPP_SRC})
add_library(class::libclass_cpp ALIAS libclass_cpp)
set_target_properties(libclass_cpp PROPERTIES CXX_STANDARD 14)
target_link_libraries(libclass_cpp PUBLIC libclass)
target_include_directories(libclass_cpp PUBLIC cpp)
set_target_properties(libclass_cpp PROPERTIES OUTPUT_NAME class_cpp)



## Python (classy)
## TODO ##
