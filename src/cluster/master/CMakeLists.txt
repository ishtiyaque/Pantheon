# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT license.

cmake_minimum_required(VERSION 3.10)

project(Pantheon VERSION 1.0 LANGUAGES CXX)

include(CheckCXXCompilerFlag)
function(seal_enable_cxx_compiler_flag_if_supported flag)
    string(FIND "${CMAKE_CXX_FLAGS}" "${flag}" flag_already_set)
    if(flag_already_set EQUAL -1)
        message(STATUS "Adding CXX compiler flag: ${flag} ...")
        check_cxx_compiler_flag("${flag}" flag_supported)
        if(flag_supported)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE)
        endif()
        unset(flag_supported CACHE)
    endif()
endfunction()



#set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -Wall -O0")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

seal_enable_cxx_compiler_flag_if_supported("-O3")
# Executable will be in ./bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

#include_directories( ${CMAKE_SOURCE_DIR}/../common)


# Import Microsoft SEAL
find_package(SEAL 0.0.0 EXACT REQUIRED)

find_package(OpenMP REQUIRED)

include_directories( ${CMAKE_SOURCE_DIR}/../common)

add_executable(master master.cpp )

# Link Microsoft SEAL
target_link_libraries(master SEAL::seal OpenMP::OpenMP_CXX rpc)


