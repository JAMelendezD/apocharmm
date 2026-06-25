add_subdirectory(unittests)

if(APOCHARMM_ENABLE_PYTHON_API)
  add_subdirectory(pytest)
endif()
