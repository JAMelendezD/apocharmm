cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

project(apoCHARMM
  VERSION 1.0.0
  LANGUAGES C CXX CUDA
)

# Only do these if this is the main project, and not if it is included through
# add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  # Generate compile_commands.json by default for tools such as clangd.
  # Useres can disable this with: -DCMAKE_EXPORT_COMPILE_COMMANDS=OFF
  #
  # If an old build directory has an empty cached value, treat that as
  # "not configured" and reset it to the project default.
  if(DEFINED CMAKE_EXPORT_COMPILE_COMMANDS AND
      "${CMAKE_EXPORT_COMPILE_COMMANDS}" STREQUAL "")
    unset(CMAKE_EXPORT_COMPILE_COMMANDS CACHE)
  endif()
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL
    "Export compile_commands.json for language servers")

  # Specify the C++ standard
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED True)

  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
  endif()

  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${apoCHARMM_SOURCE_DIR}/bin
  )

  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    ${apoCHARMM_SOURCE_DIR}/lib
  )

  # Compiler flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3 -O0")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wall -Wextra -Wno-deprecated-gpu-targets -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_90,code=sm_90")
  set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -g3 -G -O0")
  set(CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS_RELEASE} -O3 -lineinfo")
endif()

# CUDA Toolkit
find_package(CUDAToolkit REQUIRED)
set(CMAKE_CUDA_USE_RESPONSE_FILE_FOR_INCLUDES OFF)

# Build the C++/CUDA apoCHARMM library
add_subdirectory(${apoCHARMM_SOURCE_DIR}/src)

option(APOCHARMM_ENABLE_PYTHON_API "Build apoCHARMM Python interface (C ABI)" ON)
if(APOCHARMM_ENABLE_PYTHON_API)
  add_subdirectory(${apoCHARMM_SOURCE_DIR}/capi)
endif()

# Compile the test cases
add_subdirectory(${apoCHARMM_SOURCE_DIR}/test)
