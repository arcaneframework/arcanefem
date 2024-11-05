#!/bin/bash

#======================================================================================
# Script Overview:
# This script runs benchmarks for ArcaneFEM on multiple configurations.
# It measures execution times for CPU and GPU matrix assembly methods on different mesh sizes.
# The benchmark results are saved to CSV files for further analysis.
#======================================================================================

# Directory where all benchmark results will be saved
WORKING_DIR="benchmark-output"

# Path to the executable for running tests
EXECUTABLE="$HOME/Git/toutane/arcanefem_forked/build/poisson/Poisson"

# Python script that extracts specified metrics from the ArcaneFEM JSON output
# The metrics can be specified as command line arguments
PYTHON_SCRIPT="$HOME/Git/toutane/arcanefem_forked/poisson/get_stats_from_json.py"

#--------------------------------------------------------------------------------------
# Mesh configurations
# Different mesh sizes (small, medium, large) for 2D and 3D simulations
#--------------------------------------------------------------------------------------

SIZES=("small" "medium" "large")

# 2D mesh templates and paths for each size
TEMPLATE_FILENAME_2D="$HOME/Git/toutane/arcanefem_forked/poisson/TEST_TEMPLATE_2D.xml"
MESH_2D_SMALL="$HOME/Git/toutane/arcanefem_forked/poisson/meshes/L-shape-small.msh"
MESH_2D_MEDIUM="$HOME/Git/toutane/arcanefem_forked/poisson/meshes/L-shape-medium.msh"
MESH_2D_LARGE="$HOME/Git/toutane/arcanefem_forked/poisson/meshes/L-shape-large.msh"

# 3D mesh templates and paths for each size
TEMPLATE_FILENAME_3D="$HOME/Git/toutane/arcanefem_forked/poisson/TEST_TEMPLATE_3D.xml"
MESH_3D_SMALL="$HOME/Git/toutane/arcanefem_forked/poisson/meshes/L-shape-3D-small.msh"
MESH_3D_MEDIUM="$HOME/Git/toutane/arcanefem_forked/poisson/meshes/L-shape-3D-medium.msh"
MESH_3D_LARGE="$HOME/Git/toutane/arcanefem_forked/poisson/meshes/L-shape-3D-large.msh"

#--------------------------------------------------------------------------------------
# Format types
# List of formats to test on CPU and GPU, along with their display names for reporting
#--------------------------------------------------------------------------------------

# CPU matrix formats
CPU_FORMATS=("legacy" "coo" "coo-sorting" "csr")
CPU_FORMATS_MAJ=("Legacy" "Coo" "CooSort" "Csr")

# GPU matrix formats
GPU_FORMATS=("coo-gpu" "coo-sorting-gpu" "csr-gpu" "blcsr" "nwcsr")
GPU_FORMATS_MAJ=("Coo_Gpu" "CooSort_Gpu" "Csr_Gpu" "CsrBuildLess" "CsrNodeWise")

# Number of MPI instances to test for each configuration
CPU_CORE_NUMBERS=(1)

# Cache warming passed as argument to executable
CACHE_WARMING=1

#--------------------------------------------------------------------------------------
# Prepare and run tests
# The main function loop iterates over all configurations and formats to launch tests.
#--------------------------------------------------------------------------------------

# Ensure the working directory exists and switch to it
mkdir -p "$WORKING_DIR"
cd "$WORKING_DIR" || exit 1

#--------------------------------------------------------------------------------------
# launchTestCpu
# Runs a CPU test, parses results from the JSON file, and saves relevant data to CSV
# Parameters:
# - test_file: Path to the XML file with CPU test configuration
# - instance_num: Number of MPI instances
# - res_file: Path to CSV file where results are stored
#--------------------------------------------------------------------------------------
launchTestCpu() {
  local test_file=$1
  local instance_num=$2
  local res_file=$3

  # Run CPU test with MPI and save JSON results if successful
  if mpirun -n "$instance_num" "$EXECUTABLE" "$test_file" "-A,CACHE_WARMING=${CACHE_WARMING}" > /dev/null 2>&1; then
    echo -e "OK ${test_file}"
    mv "./output/listing/time_stats.json" "./time_stats.json"

    # Extract key metrics from JSON using Python script
    python "$PYTHON_SCRIPT" "./time_stats.json" "BuildMatrix,AddAndCompute" > "brief.txt"

    # Parse execution times for each format and add them to CSV
    line=$(grep "Cell" "brief.txt" | awk '{print $2}')
    for format in "${CPU_FORMATS_MAJ[@]}"; do
      time=$(grep "AssembleBilinearOperator_${format}:" "brief.txt" | awk '{print $2}')
      line="${line},${time}"
    done
    echo "$line" >> "$res_file"
  else
    echo -e "\e[31mFAIL ${test_file}"
    echo -e "command: mpirun -n ${instance_num} ${EXECUTABLE} ${test_file}\e[0m"
  fi

  # Move log file for reference
  mv "./output/listing/logs.0" "./logs.0"
}

#--------------------------------------------------------------------------------------
# launchTestGpu
# Similar to launchTestCpu, but runs tests with GPU-acceleration enabled
# Parameters:
# - test_file: Path to the XML file with GPU test configuration
# - instance_num: Number of MPI instances
# - res_file: Path to CSV file where results are stored
#--------------------------------------------------------------------------------------
launchTestGpu() {
  local test_file=$1
  local instance_num=$2
  local res_file=$3

  if mpirun -n "$instance_num" "$EXECUTABLE" "$test_file" "-A,CACHE_WARMING=${CACHE_WARMING}" "-A,AcceleratorRuntime=cuda" > /dev/null 2>&1; then
    echo -e "OK ${test_file}"
    mv "./output/listing/time_stats.json" "./time_stats.json"

    python "$PYTHON_SCRIPT" "./time_stats.json" "BuildMatrix,AddAndCompute" > "brief.txt"

    line=""
    for format in "${GPU_FORMATS_MAJ[@]}"; do
      time=$(grep "AssembleBilinearOperator_${format}:" "brief.txt" | awk '{print $2}')
      line="${line},${time}"
    done
    sed -i "$ s/$/${line}/" "$res_file"
  else
    echo -e "\e[31mFAIL ${test_file}"
    echo -e "command: mpirun -n ${instance_num} ${EXECUTABLE} ${test_file} -A,AcceleratorRuntime=cuda\e[0m"
  fi

  mv "./output/listing/logs.0" "./logs.0"
}

#--------------------------------------------------------------------------------------
# clearTest
# Cleans up output directories between tests to avoid data carryover
#--------------------------------------------------------------------------------------
clearTest() {
  rm -rf "output" "fatal_4" # fatal_4 file is always empty, maybe an error of ArcaneFEM
}

#======================================================================================
# Main test loop
# Iterates over 2D and 3D configurations, CPU core counts, and mesh sizes.
# Runs tests for each size, dimension, and format configuration, then logs results.
#======================================================================================
for dim in {2..3}; do
  dim_dir="${dim}D"
  mkdir -p "$dim_dir"
  cd "$dim_dir" || exit 1

  for cpu_n in "${CPU_CORE_NUMBERS[@]}"; do
    cpu_dir="${cpu_n}-mpi-instance"
    mkdir -p $cpu_dir
    cd "$cpu_dir" || exit 1

    res_file="${cpu_n}-mpi-instance-results.csv"

    # Add columns name in csv file
    echo "nb-cell" > "$res_file"
    for format in "${CPU_FORMATS[@]}" "${GPU_FORMATS[@]}"; do
      sed -i "1s/$/,${format}/" "$res_file"
    done

    for size in "${SIZES[@]}"; do
      size_dir="$size"
      mkdir -p "$size_dir"
      cd "$size_dir" || exit 1

      # Resolve mesh file name
      mesh_var="MESH_${dim}D_${size}"
      mesh_var="${mesh_var^^}" # To upper case
      mesh_file=${!mesh_var}

      template_var="TEMPLATE_FILENAME_${dim}D"

      # Run CPU test
      mkdir -p "cpu"
      cd "cpu" || exit 1

      test_filename="Test.${dim}D.${cpu_n}-mpi-instance.${size}.cpu-formats.arc"
      sed "s|MESH_FILE|${mesh_file}|g" "${!template_var}" > "$test_filename" # Replace mesh filename in template

      for format in "${CPU_FORMATS[@]}"; do # Add formats in template
        sed -i "/<!-- FORMATS -->/a \\
        <${format}>true</$format>" "$test_filename"
      done

      launchTestCpu "$test_filename" "$cpu_n" "../../$res_file"
      clearTest

      cd "../" # Back to size directory

      # Run GPU test, do the same as for CPU
      mkdir -p "gpu"
      cd "gpu" || exit 1

      test_filename="Test.${dim}D.${cpu_n}-mpi-instance.${size}.gpu-formats.arc"
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


