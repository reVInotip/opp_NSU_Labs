cmake_minimum_required(VERSION 3.5)

project(jacobi_algo)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)

include_directories(../include)

set(SOURCE_EXE ../src/jacobi_algo.c  ../src/utils.c)
set(SOURCE_LIB_FUNCTION ../src/function/function.c ../src/function/node_function.c)
set(SOURCE_LIB_BUFFER ../src/buffer.c)

add_library(function STATIC ${SOURCE_LIB_FUNCTION} ${SOURCE_LIB_BUFFER})
add_library(buffer STATIC ${SOURCE_LIB_BUFFER})

add_executable(${PROJECT_NAME} ${SOURCE_EXE})

find_package(MPI REQUIRED)

target_link_libraries(${PROJECT_NAME} buffer function MPI::MPI_C)
target_compile_options(${PROJECT_NAME} PRIVATE
  -Wall -Wextra -Wpedantic -O1
)