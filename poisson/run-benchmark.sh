#!/bin/bash

#======================================================================================
# Script Overview:
# This script runs benchmarks for ArcaneFEM on multiple configurations.
# It measures execution times for CPU and GPU matrix assembly methods on different mesh sizes.
# The benchmark results are saved to TSV files for further analysis.
#======================================================================================

# Path to the executable for running tests
EXECUTABLE="$(pwd)/Poisson"

# Python script that extracts specified metrics from the ArcaneFEM JSON output
# The metrics can be specified as command line arguments
PYTHON_SCRIPT="$(pwd)/get_stats_from_json.py"

#--------------------------------------------------------------------------------------
# Mesh configurations
# Different mesh sizes (small, medium, large) for 2D and 3D simulations
#--------------------------------------------------------------------------------------

SIZES=("small" "medium" "large")

DIMENSIONS=(2 3)

# 2D mesh templates and paths for each size
TEMPLATE_FILENAME_2D="$(pwd)/TEST_TEMPLATE_2D.xml"
MESH_2D_SMALL="$(pwd)/circle_cut-small.msh"
MESH_2D_MEDIUM="$(pwd)/circle_cut-medium.msh"
MESH_2D_LARGE="$(pwd)/circle_cut-large.msh"

# 3D mesh templates and paths for each size
TEMPLATE_FILENAME_3D="$(pwd)/TEST_TEMPLATE_3D.xml"
MESH_3D_SMALL="$(pwd)/sphere_cut-small.msh"
MESH_3D_MEDIUM="$(pwd)/sphere_cut-medium.msh"
MESH_3D_LARGE="$(pwd)/sphere_cut-large.msh"

#--------------------------------------------------------------------------------------
# Format types
# List of formats to test on CPU and GPU, along with their display names for reporting
#--------------------------------------------------------------------------------------

# Formats to test
# Attention ! When using Hypre linear system, "legacy" format won't work
# blcsr is the last used format in ArcaneFEM
CPU_FORMATS=("legacy" "coo" "csr")
GPU_FORMATS=("coo-gpu" "csr-gpu" "nwcsr" "blcsr")

# Number of MPI instances to test for each configuration
CPU_CORE_NUMBERS=(1 2)

# Cache warming passed as argument to executable
CACHE_WARMING=10

# Directory where all benchmark results will be saved
DATE=$(date '+%Y-%m-%d_%H:%M:%S')
WORKING_DIR="benchmark-output_${CACHE_WARMING}-cw_${DATE}"

# Arcane AcceleratorRuntime passed as argument to executable
ACCELERATOR_RUNTIME="cuda"

#--------------------------------------------------------------------------------------
# Prepare and run tests
# The main function loop iterates over all configurations and formats to launch tests.
#--------------------------------------------------------------------------------------

# Ensure the working directory exists and switch to it
mkdir -p "$WORKING_DIR"
cd "$WORKING_DIR" || exit 1

#--------------------------------------------------------------------------------------
# Runs a CPU test for CPU and GPU formats (with GPU-acceleration disabled), parses results
# from the JSON file, and saves relevant data to TSV
# Parameters:
# - test_file: Path to the XML file with CPU test configuration
# - instance_num: Number of MPI instances
# - res_file: Path to TSV file where results are stored
#--------------------------------------------------------------------------------------
launchTestCpu() {
  local test_file=$1
  local instance_num=$2
  local res_file=$3

  if [ ! -e "$EXECUTABLE" ]; then
    echo -e "\e[31mExecutable file: \"${EXECUTABLE}\" not found, stop\e[0m"
    exit
  fi

  if [ ! -e "$test_file" ]; then
    echo -e "\e[31mTest file: \"${test_file}\" not found, stop\e[0m"
    exit
  fi

  # Run CPU test with MPI and save JSON results if successful
  echo "Info: Starting ${test_file}..."
  if mpirun -n "$instance_num" "$EXECUTABLE" "$test_file" "-A,CACHE_WARMING=${CACHE_WARMING}" > "stdout.txt" 2> "stderr.txt"; then
    mv "./output/listing/time_stats.json" "./time_stats.json"

    # Extract key metrics from JSON using Python script
    if python "$PYTHON_SCRIPT" "./time_stats.json" "BuildMatrix,AddAndCompute" > "brief.txt"; then

      # Parse execution times for each format and add them to TSV
      line=$(grep "Element" "brief.txt" | awk '{print $2}')

      for format in "${ALL_CPU_FORMATS[@]}"; do
        if contains "$format" "${CPU_FORMATS[@]}"; then
          time=$(grep "AssembleBilinearOperator_${format}:" "brief.txt" | awk '{print $2}')
          line+="\t${time}"
        else
          line+="\tNaN"
        fi
      done

      for format in "${ALL_GPU_FORMATS[@]}"; do
        if contains "$format" "${GPU_FORMATS[@]}"; then
          time=$(grep "AssembleBilinearOperator_${format}:" "brief.txt" | awk '{print $2}')
          line+="\t${time}"
        else
          line+="\tNaN"
        fi
      done

      echo -e "$line" >> "$res_file"

      mv "./output/listing/logs.0" "./logs.0"
      echo -e "Info: Done\n"
    else
      echo -e "\e[31mAn error occured in ${PYTHON_SCRIPT}, stop\e[0m"
      exit
    fi
  else
    echo -e "\e[31mFAIL ${test_file} (command was: mpirun -n ${instance_num} ${EXECUTABLE} ${test_file} -A,CACHE_WARMING=${CACHE_WARMING}), stop\e[0m"
    exit
  fi
}

#--------------------------------------------------------------------------------------
# Similar to launchTestCpu, but runs tests for GPU format with GPU-acceleration enabled
# Parameters:
# - test_file: Path to the XML file with GPU test configuration
# - instance_num: Number of MPI instances
# - res_file: Path to TSV file where results are stored
#--------------------------------------------------------------------------------------
launchTestGpu() {
  local test_file=$1
  local instance_num=$2
  local res_file=$3

  if [ ! -e "$EXECUTABLE" ]; then
    echo -e "\e[31mExecutable file: \"${EXECUTABLE}\" not found, stop\e[0m"
    exit
  fi

  if [ ! -e "$test_file" ]; then
    echo -e "\e[31mTest file: \"${test_file}\" not found, stop\e[0m"
    exit
  fi

  echo "Info: Starting ${test_file}..."
  if mpirun -n "$instance_num" "$EXECUTABLE" "$test_file" "-A,CACHE_WARMING=${CACHE_WARMING}" "-A,AcceleratorRuntime=${ACCELERATOR_RUNTIME}" > "stdout.txt" 2> "stderr.txt"; then
    mv "./output/listing/time_stats.json" "./time_stats.json"

    if python "$PYTHON_SCRIPT" "./time_stats.json" "BuildMatrix,AddAndCompute" > "brief.txt"; then

      line=""

      for format in "${ALL_GPU_FORMATS[@]}"; do
        if contains "$format" "${GPU_FORMATS[@]}"; then
          time=$(grep "AssembleBilinearOperator_${format}:" "brief.txt" | awk '{print $2}')
          line+="\t${time}"
        else
          line+="\tNaN"
        fi
      done

      sed -i "$ s/$/${line}/" "$res_file"
      mv "./output/listing/logs.0" "./logs.0"
      echo -e "Info: Done\n"
    else
      echo -e "\e[31mAn error occured in ${PYTHON_SCRIPT}, stop\e[0m"
      exit
    fi
  else
    echo -e "\e[31mFAIL ${test_file} (command was: mpirun -n ${instance_num} ${EXECUTABLE} ${test_file} -A,CACHE_WARMING=${CACHE_WARMING} -A,AcceleratorRuntime=${ACCELERATOR_RUNTIME}), stop\e[0m"
    exit
  fi
}

#--------------------------------------------------------------------------------------
# Cleans up output directories between tests to avoid data carryover
#--------------------------------------------------------------------------------------
clearTest() {
  rm -rf "output" "fatal_4" # fatal_4 file is always empty, maybe an error of ArcaneFEM
}

#--------------------------------------------------------------------------------------
# 0：match, 1: failed
#--------------------------------------------------------------------------------------
contains() {
  local element="$1"
  shift
  for i in "$@"; do
    if [[ $i == "$element" ]]; then
      return 0
    fi
  done
  return 1
}

#--------------------------------------------------------------------------------------
# All formats supported by ArcaneFEM (do not modify)
#--------------------------------------------------------------------------------------
ALL_CPU_FORMATS=("legacy" "coo" "coo-sorting" "csr")
ALL_GPU_FORMATS=("coo-gpu" "coo-sorting-gpu" "csr-gpu" "blcsr" "nwcsr")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

function printInfo() {
  echo -n "Info: Formats: "
  for f in "${CPU_FORMATS[@]}" "${GPU_FORMATS[@]}"; do
    echo -n "| ${f} "
  done
  echo "|"
  echo "Info: Cache warming: ${CACHE_WARMING}"
  echo "Info: Output directory: '${WORKING_DIR}'"
  echo "Info: CT (cpu-thread) and CTG (cpu-thread-gpu) formats will be run without accelerator"
  echo "Info: CTG formats will be then run with AcceleratorRuntime: '${ACCELERATOR_RUNTIME}'"
  echo "-----------------------------------------------------------------------------------------"
}

#======================================================================================
# Main test loop
# Iterates over 2D configurations, CPU core counts, and mesh sizes.
# Runs tests for each size, dimension, and format configuration, then logs results.
#======================================================================================

printInfo

for dim in "${DIMENSIONS[@]}"; do
  echo -e "Info: Starting ${dim}D meshes...\n"

  dim_dir="${dim}D"
  mkdir -p "$dim_dir"
  cd "$dim_dir" || exit 1

  for cpu_n in "${CPU_CORE_NUMBERS[@]}"; do
    cpu_dir="${cpu_n}-mpi-instance"
    mkdir -p $cpu_dir
    cd "$cpu_dir" || exit 1

    res_file="results.tsv"

    # Add columns name in tsv file
    output="#nb-elt"

    for format in "${ALL_CPU_FORMATS[@]}"; do
      output+="\t${format}"
    done

    for format in "${ALL_GPU_FORMATS[@]}"; do
      output+="\t${format}-cpu"
    done

    for format in "${ALL_GPU_FORMATS[@]}"; do
      output+="\t${format}"
    done

    echo -e "$output" > "$res_file"

    for size in "${SIZES[@]}"; do
      size_dir="$size"
      mkdir -p "$size_dir"
      cd "$size_dir" || exit 1

      # Resolve mesh file name
      mesh_var="MESH_${dim}D_${size}"
      mesh_var="${mesh_var^^}" # To upper case
      mesh_file=${!mesh_var}

      if [ ! -e "$mesh_file" ]; then
        echo -e "\e[31mMeshfile: \"${mesh_file}\" not found, stop\e[0m"
        exit
      fi

      template_var="TEMPLATE_FILENAME_${dim}D"

      # Run CPU test
      mkdir -p "cpu"
      cd "cpu" || exit 1

      test_filename="Test.${dim}D.${cpu_n}-mpi-instance.${size}.accelerator-disable.arc"
      sed "s|MESH_FILE|${mesh_file}|g" "${!template_var}" > "$test_filename" # Replace mesh filename in template

      for format in "${CPU_FORMATS[@]}" "${GPU_FORMATS[@]}"; do # Add formats in template
        sed -i "/<!-- FORMATS -->/a \\
        <${format}>true</$format>" "$test_filename"
      done

      launchTestCpu "$test_filename" "$cpu_n" "../../$res_file"
      clearTest

      cd "../" # Back to size directory

      # Run GPU test, do the same as for CPU
      mkdir -p "gpu"
      cd "gpu" || exit 1

      test_filename="Test.${dim}D.${cpu_n}-mpi-instance.${size}.accelerator-enable.arc"
      sed "s|MESH_FILE|${mesh_file}|g" "${!template_var}" > "$test_filename"

      for format in "${GPU_FORMATS[@]}"; do
        sed -i "/<!-- FORMATS -->/a \\
        <${format}>true</$format>" "$test_filename"
      done

      launchTestGpu "$test_filename" "$cpu_n" "../../$res_file"
      clearTest

      cd "../" # Back to size directory
      cd "../" # Back to mpi-instance directory
    done
    cd "../" # Back to dimension directory
  done
  cd "../" # Back to working directory
done
