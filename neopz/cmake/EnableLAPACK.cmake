function(enable_lapack target)
  if(NOT USING_MKL)
    if(NOT APPLE) # In apple computer, we then use Accelerate framework
      find_package(LAPACK REQUIRED)
      target_link_libraries(${target} PRIVATE ${LAPACK_LIBRARIES} )
      foreach(LA_LIB ${LAPACK_LIBRARIES})
        if(LA_LIB MATCHES ".*mkl.*")
          target_compile_definitions(${target} PRIVATE MKLLAPACK)
          get_filename_component(LA_LIB_DIR ${LA_LIB} DIRECTORY)
          if(MKL_INCLUDE_DIR)
            set(LAPACK_INCLUDE_DIR ${MKL_INCLUDE_DIR})
            break()
          endif()
          include(cmake/FindIncludeFromMKL.cmake)
          find_include_from_mkl(LAPACK_INCLUDE_DIR ${LA_LIB_DIR} mkl_lapacke.h)
          if(NOT LAPACK_INCLUDE_DIR)
            message(FATAL_ERROR "Could not find LAPACK include directory")
          else()
            message(STATUS "LAPACK include dirs: ${LAPACK_INCLUDE_DIR}")
            target_include_directories(${target} PRIVATE ${LAPACK_INCLUDE_DIR})
            mark_as_advanced(LAPACK_INCLUDE_DIR)
          endif()
          break()
        endif()
      endforeach()
    endif()
  else()
    target_compile_definitions(${target} PRIVATE MKLLAPACK)
  endif()
  target_compile_definitions(${target} PRIVATE USING_LAPACK)
  target_compile_definitions(${target} INTERFACE PZ_USING_LAPACK)
endfunction()
