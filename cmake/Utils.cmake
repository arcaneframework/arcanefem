
unset(ARCANE_ACCELERATOR_RUNTIME)
if(ARCANE_HAS_ACCELERATOR)
  if(ARCANE_HAS_CUDA)
    message(STATUS "Arcane has been compiled with CUDA.")
    set(ARCANE_ACCELERATOR_RUNTIME "cuda")
  endif()
  if(ARCANE_HAS_HIP)
    message(STATUS "Arcane has been compiled with ROCM/HIP.")
    set(ARCANE_ACCELERATOR_RUNTIME "hip")
  endif()
endif()

# ----------------------------------------------------------------------------
# Usage:
#
# arcanefem_add_gpu_test(NAME test_name [NB_MPI n] COMMAND exe_fileARGS exe_args)
#
macro(arcanefem_add_gpu_test)
  set(options)
  set(oneValueArgs NB_MPI NAME COMMAND)
  set(multiValueArgs ARGS)
  cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if (NOT ARGS_NAME)
    message(FATAL_ERROR "No arg NAME for macro 'arcanefem_add_gpu_test'")
  endif()
  if (NOT ARGS_COMMAND)
    message(FATAL_ERROR "No arg COMMAND for macro 'arcanefem_add_gpu_test'")
  endif()
  # Add test without accelerator Runtime
  if (ARGS_NB_MPI)
    if (MPIEXEC_EXECUTABLE)
      add_test(NAME ${ARGS_NAME} COMMAND ${MPIEXEC_EXECUTABLE} -n ${ARGS_NB_MPI} ${ARGS_COMMAND} ${ARGS_ARGS})
    endif()
  else()
    add_test(NAME ${ARGS_NAME} COMMAND ${ARGS_COMMAND} ${ARGS_ARGS})
  endif()

  # Add test WITH accelerator Runtime
  if(ARCANE_HAS_ACCELERATOR)
    set(_RUNTIME_ARGS "-A,AcceleratorRuntime=${ARCANE_ACCELERATOR_RUNTIME}")
    if (ARGS_NB_MPI)
      if (MPIEXEC_EXECUTABLE)
        add_test(NAME ${ARGS_NAME}_${ARCANE_ACCELERATOR_RUNTIME} COMMAND ${MPIEXEC_EXECUTABLE} -n ${ARGS_NB_MPI} ${ARGS_COMMAND} ${_RUNTIME_ARGS} ${ARGS_ARGS})
      endif()
    else()
      add_test(NAME ${ARGS_NAME}_${ARCANE_ACCELERATOR_RUNTIME} COMMAND ${ARGS_COMMAND} ${_RUNTIME_ARGS} ${ARGS_ARGS})
    endif()
  endif()
endmacro()

# ----------------------------------------------------------------------------
# Local Variables:
# tab-width: 2
# indent-tabs-mode: nil
# coding: utf-8-with-signature
# End:
