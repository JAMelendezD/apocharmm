# Add any new tests to this list!
set(SOURCES_UNITTESTS_TESTS
  # Base tests
  unittest-deviceVector.cpp
  unittest-cudaContainer.cpp

  # # Manager tests
  # unittest-forceManager.cpp

  # Object tests
  unittest-charmmParameters.cpp
  unittest-charmmPSF.cpp
  unittest-coordinates.cpp
  unittest-charmmCrd.cpp
  unittest-atomSelection.cpp

  # Subscriber tests
  unittest-dcdSubscriber.cpp
  unittest-restartSubscriber.cpp

  # Integrator tests
  unittest-cudaIntegrator.cpp
  unittest-cudaNoseHooverIntegrator.cpp
  unittest-cudaLangevinThermostatIntegrator.cpp
  unittest-cudaLangevinPistonIntegrator.cpp

  # Restraint tests
  unittest-harmonicRestraintForce.cpp
  unittest-harmonicCenterOfMassRestraintForce.cpp
)

#
# Stuff below shouldn't need to be modified
#
set(SOURCES_UNITTESTS_MAIN unittest-main.cpp)

# A function to check that the given hard-coded list contains all headers + sources in the given folder.
function(CheckFileList LIST_VAR FOLDER)
  set(MESSAGE " should be added to the variable ${LIST_VAR}")
  set(MESSAGE "${MESSAGE} in ${CMAKE_CURRENT_LIST_FILE}\n")
  file(GLOB GLOBBED_LIST "${FOLDER}/*.cpp"
    "${FOLDER}/*.hpp")
  list(REMOVE_ITEM GLOBBED_LIST ${${LIST_VAR}})
endfunction()

set(UNITTEST_DIR ${PROJECT_SOURCE_DIR}/test/unittests)

set(HEADER_DIR ${UNITTEST_DIR})

# Check if all sources are listed; warn if not.
set(SOURCES_ALL
  ${SOURCES_UNITTESTS_MAIN}
  ${SOURCES_UNITTESTS_TESTS}
)

foreach(name ${SOURCES_ALL})
    list(APPEND SOURCES_ALL_PATH ${UNITTEST_DIR}/${name})
endforeach()

CheckFileList(SOURCES_ALL_PATH ${UNITTEST_DIR})

# Create serial target names.
string(REPLACE ".cpp" "" BASENAMES_UNITTESTS_TESTS "${SOURCES_UNITTESTS_TESTS}")
string(REPLACE ".cu" "" BASENAMES_UNITTESTS_TESTS "${BASENAMES_UNITTESTS_TESTS}")
set(TARGETS_UNITTESTS_TESTS ${BASENAMES_UNITTESTS_TESTS})

set(TARGETS_ALL ${TARGETS_UNITTESTS_TESTS} CatchMain)

# Define program targets:
add_library(CatchMain OBJECT ${UNITTEST_DIR}/${SOURCES_UNITTESTS_MAIN} ${HEADER_DIR}/catch.hpp)

foreach(name ${TARGETS_UNITTESTS_TESTS})
  set(filename ${name}.cpp)

  add_executable(${name} ${filename} $<TARGET_OBJECTS:CatchMain> ${HEADER_DIR}/catch.hpp)

  target_compile_definitions(${name} PUBLIC DATA_PATH="${PROJECT_SOURCE_DIR}/test/data/")
  target_link_libraries(${name} PRIVATE
    ${EXTERNAL_LIBRARIES}
    apoCHARMMlib
    CUDA::cudart
  )
  set_property(TARGET ${name} PROPERTY CUDA_RUNTIME_LIBRARY SHARED) # For CMake 3.23 and later

  set_target_properties(${name} PROPERTIES LINKER_LANGUAGE CUDA)
  add_test(NAME ${name}
    COMMAND ${name}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test/unittests)
endforeach()

foreach(name ${TARGETS_ALL})
  target_include_directories(${name} PRIVATE ${HEADER_DIR})

  # GCC specific warning go here
  if(CMAKE_CXX_COMPILER_ID MATCHES "GCC")
    # Actually keep these
    target_compile_options(${name} PRIVATE -Wweak-vtables -Wexit-time-destructors -Wglobal-constructors -Wmissing-noreturn -Wno-int-in-bool-context)
  endif()

  # Clang specific warning go here
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # Actually keep these
    target_compile_options(${name} PRIVATE -Wweak-vtables -Wexit-time-destructors -Wglobal-constructors)
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(${name} PRIVATE /W4 /w44265 /WX)
  endif()
endforeach()
