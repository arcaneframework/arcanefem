#!/bin/bash

# Constants
WORKING_DIR="benchmark"
EXECUTABLE="$HOME/Git/toutane/arcanefem_forked/build/poisson/Poisson"
PYTHON_SCRIPT="$HOME/Git/toutane/arcanefem_forked/poisson/get_stats_from_json.py"
TEMPLATE_FILENAME="$HOME/Git/toutane/arcanefem_forked/poisson/TEST_TEMPLATE.xml"

# Mesh files (update paths as necessary)
MESH_2D_SMALL="$HOME/Git/toutane/arcanefem_forked/meshes/msh/L-shape.msh"
MESHES=("MESH_2D_SMALL") # Add other mesh variables if needed

CPU_FORMATS=("legacy" "coo" "coo-sorting" "csr")
CPU_FORMATS_MAJ=("Legacy" "Coo" "CooSort" "Csr")

GPU_FORMATS=("coo-gpu" "coo-sorting-gpu" "csr-gpu" "blcsr" "nwcsr")
GPU_FORMATS_MAJ=("Coo_Gpu" "CooSort_Gpu" "Csr_Gpu" "CsrBuildLess" "CsrNodeWise")

SIZES=("small")
CPU_CORE_NUMBERS=(1)

# Create working directory if it doesn't exist
mkdir -p "$WORKING_DIR"
cd "$WORKING_DIR" || exit 1

# Launch test for CPU
launchTestCpu() {
  local test_file=$1
  local instance_num=$2
  local res_file=$3

  # Run test
  if mpirun -n "$instance_num" "$EXECUTABLE" "$test_file" > /dev/null 2>&1; then
    echo -e "${test_file} [${instance_num} CPU]:       OK"
    cp "./output/listing/time_stats.json" "./cpu_formats_time_stats.json"

    python "$PYTHON_SCRIPT" "./output/listing/time_stats.json" "BuildMatrix,AddAndCompute" > "cpu_formats_brief.txt"

    grep "Cell" "cpu_formats_brief.txt" | awk '{print $2}' >> $res_file

    for format in "${CPU_FORMATS_MAJ[@]}"; do
      time=$(grep "AssembleBilinearOperator_${format}:" "cpu_formats_brief.txt" | awk '{print $2}')
      sed -i "2s/$/,${time}/" "$res_file"
    done

  else
    echo -e "\e[31m${test_file} [${instance_num} CPU]:       FAIL"
    echo -e "command: mpirun -n ${instance_num} ${EXECUTABLE} ${test_file}\e[0m"
  fi

  cp "./output/listing/logs.0" "./cpu_formats_logs.0"
}

# Launch test for GPU
launchTestGpu() {
  local test_file=$1
  local instance_num=$2
  local res_file=$3

  # Run test
  if mpirun -n "$instance_num" "$EXECUTABLE" "$test_file" "-A,AcceleratorRuntime=cuda" > /dev/null 2>&1; then
    echo -e "${test_file} [${instance_num} CPU + GPU]: OK"
    cp "./output/listing/time_stats.json" "./gpu_formats_time_stats.json"

    python "$PYTHON_SCRIPT" "./output/listing/time_stats.json" "BuildMatrix,AddAndCompute" > "gpu_formats_brief.txt"

    for format in "${GPU_FORMATS_MAJ[@]}"; do
      time=$(grep "AssembleBilinearOperator_${format}:" "gpu_formats_brief.txt" | awk '{print $2}')
      sed -i "2s/$/,${time}/" "$res_file"
    done

  else
    echo -e "\e[31m${test_file} [${instance_num} CPU + GPU]: FAIL"
    echo -e "command: mpirun -n ${instance_num} ${EXECUTABLE} ${test_file} -A,AcceleratorRuntime=cuda\e[0m"
  fi

  cp "./output/listing/logs.0" "./gpu_formats_logs.0"
}

# Clear previous output
clearTest() {
  rm -rf "output" "fatal_4"
}

# Main test loop
for dim in 2; do
  dim_dir="${dim}D"
  mkdir -p "$dim_dir"
  cd "$dim_dir" || exit 1

  for cpu_n in "${CPU_CORE_NUMBERS[@]}"; do
    cpu_dir="${cpu_n}-mpi-instance"
    mkdir -p $cpu_dir
    cd "$cpu_dir" || exit 1

    # Create CSV result file
    res_file="${cpu_n}-mpi-instance-results.csv"
    echo "nb-cell" > "$res_file"

    for size in "${SIZES[@]}"; do
      size_dir="$size"
      mkdir -p "$size_dir"
      cd "$size_dir" || exit 1

      dir="cpu"
      mkdir -p "$dir"
      cd "$dir" || exit 1

      # Prepare CPU test file
      test_name="Test.${dim}D.${size}.cpu_${cpu_n}"
      mesh_file="${!MESHES[0]}" # Expand mesh variable dynamically
      test_filename="${test_name}.cpu_formats.arc"

      sed "s|MESH_FILE|${mesh_file}|g" "$TEMPLATE_FILENAME" > "$test_filename"

      for format in "${CPU_FORMATS[@]}"; do
        sed -i "1s/$/,${format}/" "../../$res_file"

        sed -i "/<!-- FORMATS -->/a \\
        <${format}>true</$format>" "$test_filename"
      done

      launchTestCpu "$test_filename" "$cpu_n" "../../$res_file"
      clearTest
      cd "../" # Back to size directory

      # Prepare GPU test file
      dir="gpu"
      mkdir -p "$dir"
      cd "$dir" || exit 1

      test_filename="${test_name}.gpu_formats.arc"
      sed "s|MESH_FILE|${mesh_file}|g" "$TEMPLATE_FILENAME" > "$test_filename"

      for format in "${GPU_FORMATS[@]}"; do
        sed -i "1s/$/,${format}/" "../../$res_file"

        sed -i "/<!-- FORMATS -->/a \\
        <${format}>true</$format>" "$test_filename"
      done

      launchTestGpu "$test_filename" "$cpu_n" "../../$res_file"
      clearTest
      cd "../.." # Back to dimension directory
    done
    # echo '' # new line
    cd "../" # Back to working directory
  done
done

