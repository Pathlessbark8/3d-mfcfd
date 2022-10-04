#MFCFD CMAKE SCRIPT

cmake_minimum_required(VERSION 3.10)
project(mfcfd)
ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if( "${MFCFD}" MATCHES "CUDA" )
ENABLE_LANGUAGE(CUDA)
include_directories(${CUDA_INCLUDE_DIRS})
set(LIBS ${LIBS} ${CUDA_LIBRARIES})

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 80)
endif()

add_compile_options(-O3 --generate-line-info -lineinfo -mcmodel=medium )
list(APPEND MFCFD_CUDA_SOURCES ${CMAKE_SOURCE_DIR}/test.cu)
set(CMAKE_C_FLAGS "${CMAKE_CUDA_FLAGS} -Mcuda=nordc,lineinfo  -Minfo")
set(CUDA_SEPARABLE_COMPILATION ON)
add_executable(execname ${MFCFD_CUDA_SOURCES})
target_link_libraries(execname PUBLIC ${LIBS})
endif()