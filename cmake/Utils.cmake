
unset(ARCANE_ACCELERATOR_RUNTIME)
if (ARCANE_HAS_ACCELERATOR)
  if (ARCANE_HAS_CUDA)
    message(STATUS "Arcane has been compiled with CUDA.")
    set(ARCANE_ACCELERATOR_RUNTIME "cuda")
  endif ()
  if (ARCANE_HAS_HIP)
    message(STATUS "Arcane has been compiled with ROCM/HIP.")
    set(ARCANE_ACCELERATOR_RUNTIME "hip")
  endif ()
endif ()

# ----------------------------------------------------------------------------
# Usage:
#
# arcanefem_add_seq_test(GPU NAME test_name MODULE module_name COMMAND exe_file ARGS exe_args)
# ("GPU" if you want GPU tests)
#
macro(arcanefem_add_seq_test)
  set(options GPU)
  set(oneValueArgs NAME COMMAND MODULE)
  set(multiValueArgs ARGS)
  cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if (NOT ARGS_NAME)
    message(FATAL_ERROR "No arg NAME for macro 'arcanefem_add_seq_test'")
  endif ()
  if (NOT ARGS_COMMAND)
    message(FATAL_ERROR "No arg COMMAND for macro 'arcanefem_add_seq_test'")
  endif ()
  if (NOT ARGS_MODULE)
    if (NOT MODULE_NAME)
      message(FATAL_ERROR "No arg MODULE or var MODULE_NAME for macro 'arcanefem_add_seq_test'")
    else ()
      set(ARGS_MODULE ${MODULE_NAME})
    endif ()
  endif ()
  if(NOT ARGS_ARGS)
    message(FATAL_ERROR "No arg ARGS for macro 'arcanefem_add_seq_test'")
  endif ()

  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_outputs/${ARGS_MODULE}")

  set(_RUNTIME_ARGS "-A,OutputDirectory=${CMAKE_CURRENT_BINARY_DIR}/test_outputs/${ARGS_MODULE}/${ARGS_NAME}")

    # --- Serial CPU test (always) ---
  #message(STATUS "Add sequential test [${ARGS_MODULE}]${ARGS_NAME}")
  add_test(NAME "[${ARGS_MODULE}]${ARGS_NAME}"
    COMMAND ${ARGS_COMMAND} ${_RUNTIME_ARGS} ${ARGS_ARGS})

  # --- Serial GPU test (if ARGS_GPU and if accelerator available) ---
  if (ARGS_GPU AND ARCANE_HAS_ACCELERATOR)
    set(_RUNTIME_ARGS "-A,OutputDirectory=${CMAKE_CURRENT_BINARY_DIR}/test_outputs/${ARGS_MODULE}/${ARGS_NAME}_${ARCANE_ACCELERATOR_RUNTIME},AcceleratorRuntime=${ARCANE_ACCELERATOR_RUNTIME}")

    #message(STATUS "Add sequential GPU test [${ARGS_MODULE}]${ARGS_NAME}_${ARCANE_ACCELERATOR_RUNTIME}")
    add_test(NAME "[${ARGS_MODULE}]${ARGS_NAME}_${ARCANE_ACCELERATOR_RUNTIME}"
      COMMAND ${ARGS_COMMAND} ${_RUNTIME_ARGS} ${ARGS_ARGS})
  endif ()
endmacro()

# ----------------------------------------------------------------------------
# Usage:
#
# arcanefem_add_mpi_test(GPU NAME test_name MODULE module_name NB_MPI test_nb_proc COMMAND exe_file ARGS exe_args)
# ("GPU" if you want GPU tests)
#
macro(arcanefem_add_mpi_test)
  set(options GPU)
  set(oneValueArgs NB_MPI NAME COMMAND MODULE)
  set(multiValueArgs ARGS)
  cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if (NOT ARGS_NAME)
    message(FATAL_ERROR "No arg NAME for macro 'arcanefem_add_mpi_test'")
  endif ()
  if (NOT ARGS_COMMAND)
    message(FATAL_ERROR "No arg COMMAND for macro 'arcanefem_add_mpi_test'")
  endif ()
  if (NOT ARGS_MODULE)
    if (NOT MODULE_NAME)
      message(FATAL_ERROR "No arg MODULE or var MODULE_NAME for macro 'arcanefem_add_mpi_test'")
    else ()
      set(ARGS_MODULE ${MODULE_NAME})
    endif ()
  endif ()
  if(NOT ARGS_ARGS)
    message(FATAL_ERROR "No arg ARGS for macro 'arcanefem_add_seq_test'")
  endif ()

  if (MPIEXEC_EXECUTABLE)

    # Default number of MPI procs for the parallel variant if not specified
    if (NOT ARGS_NB_MPI)
      set(ARGS_NB_MPI 2)
    endif ()

    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_outputs/${ARGS_MODULE}")

    # --- Parallel CPU test (always, if mpiexec available) ---
    set(_RUNTIME_ARGS "-A,OutputDirectory=${CMAKE_CURRENT_BINARY_DIR}/test_outputs/${ARGS_MODULE}/${ARGS_NAME}_${ARGS_NB_MPI}p")

    #message(STATUS "Add mpi test [${ARGS_MODULE}]${ARGS_NAME}_${ARGS_NB_MPI}p")
    add_test(NAME "[${ARGS_MODULE}]${ARGS_NAME}_${ARGS_NB_MPI}p"
      COMMAND ${MPIEXEC_EXECUTABLE} -n ${ARGS_NB_MPI} ${ARGS_COMMAND} ${_RUNTIME_ARGS} ${ARGS_ARGS})

    # --- Parallel GPU test (if mpiexec AND if ARGS_GPU and accelerator available) ---
    if (ARGS_GPU AND ARCANE_HAS_ACCELERATOR)

      set(_RUNTIME_ARGS "-A,OutputDirectory=${CMAKE_CURRENT_BINARY_DIR}/test_outputs/${ARGS_MODULE}/${ARGS_NAME}_${ARGS_NB_MPI}p_${ARCANE_ACCELERATOR_RUNTIME},AcceleratorRuntime=${ARCANE_ACCELERATOR_RUNTIME}")

      #message(STATUS "Add mpi GPU test [${ARGS_MODULE}]${ARGS_NAME}_${ARGS_NB_MPI}p_${ARCANE_ACCELERATOR_RUNTIME}")
      add_test(NAME "[${ARGS_MODULE}]${ARGS_NAME}_${ARGS_NB_MPI}p_${ARCANE_ACCELERATOR_RUNTIME}"
        COMMAND ${MPIEXEC_EXECUTABLE} -n ${ARGS_NB_MPI} ${ARGS_COMMAND} ${_RUNTIME_ARGS} ${ARGS_ARGS})
    endif ()
  endif ()
endmacro()

# ----------------------------------------------------------------------------
# Usage:
#
# arcanefem_add_full_test(NAME test_name MODULE module_name NB_MPI test_nb_proc COMMAND exe_file ARGS exe_args)
#
macro(arcanefem_add_full_test)
  set(options)
  set(oneValueArgs NB_MPI NAME COMMAND MODULE)
  set(multiValueArgs ARGS)
  cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if (NOT ARGS_NAME)
    message(FATAL_ERROR "No arg NAME for macro 'arcanefem_add_full_test'")
  endif ()
  if (NOT ARGS_COMMAND)
    message(FATAL_ERROR "No arg COMMAND for macro 'arcanefem_add_full_test'")
  endif ()
  if (NOT ARGS_MODULE)
    if (NOT MODULE_NAME)
      message(FATAL_ERROR "No arg MODULE or var MODULE_NAME for macro 'arcanefem_add_full_test'")
    else ()
      set(ARGS_MODULE ${MODULE_NAME})
    endif ()
  endif ()
  if(NOT ARGS_ARGS)
    message(FATAL_ERROR "No arg ARGS for macro 'arcanefem_add_seq_test'")
  endif ()

  # Default number of MPI procs for the parallel variant if not specified
  if (NOT ARGS_NB_MPI)
    set(ARGS_NB_MPI 2)
  endif ()

  arcanefem_add_seq_test(GPU NAME ${ARGS_NAME}
    MODULE ${ARGS_MODULE}
    COMMAND ${ARGS_COMMAND}
    ARGS ${ARGS_ARGS})

  arcanefem_add_mpi_test(GPU NAME ${ARGS_NAME}
    MODULE ${ARGS_MODULE}
    NB_MPI ${ARGS_NB_MPI}
    COMMAND ${ARGS_COMMAND}
    ARGS ${ARGS_ARGS})

endmacro()

# ----------------------------------------------------------------------------
# Local Variables:
# tab-width: 2
# indent-tabs-mode: nil
# coding: utf-8-with-signature
# End:
