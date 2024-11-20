#!/bin/bash

#======================================================================================
# This script runs benchmarks for ArcaneFEM on multiple configurations.
# It measures execution times for CPU and GPU bilinear assembly methods on different
# mesh sizes. The benchmark results are saved to TSV files for further analysis.
#======================================================================================

EXECUTABLE="$(pwd)/Poisson"
PYTHON_SCRIPT="$(pwd)/get_stats_from_json.py"

CACHE_WARMING=5
WORKING_DIR_BASE="benchmark-output.${CACHE_WARMING}-cw"
ACCELERATOR_RUNTIME="cuda"

MPI_N=(1 2 4 8)
CPU_FORMATS=("legacy" "coo" "csr")

MPI_N_ACCELERATED=(1 2 4 8)
GPU_FORMATS=("coo-gpu" "csr-gpu" "nwcsr" "blcsr")

DIMENSIONS=(2 3)
TEMPLATE_FILENAMES=("$(pwd)/TEST_TEMPLATE_2D.xml" "$(pwd)/TEST_TEMPLATE_3D.xml")
SIZES=("small" "medium" "large")

declare -A MESHES_2D=(
["small"]="$(pwd)/circle_cut-small.msh"
["medium"]="$(pwd)/circle_cut-medium.msh"
["large"]="$(pwd)/circle_cut-large.msh"
)

declare -A MESHES_3D=(
["small"]="$(pwd)/sphere_cut-small.msh"
["medium"]="$(pwd)/sphere_cut-medium.msh"
["large"]="$(pwd)/sphere_cut-large.msh"
)

# Adastra configuration
IS_ADASTRA=false # replace by `true` to apply
$IS_ADASTRA && ACCELERATOR_RUNTIME="hip"
CPU_PER_TASK=1 # no effect if -T option is disabled
THREADS_PER_CORE=1 # disable hyper-threading
SRUN_ARGS="--cpus-per-task=${CPU_PER_TASK} --threads-per-core=${THREADS_PER_CORE}"
SRUN_ARGS_TO_ADD_IF_ACCELERATED="--gpus=1"

#======================================================================================
# Utility Functions
#======================================================================================

ALL_CPU_FORMATS=("legacy" "coo" "coo-sorting" "csr")
ALL_GPU_FORMATS=("coo-gpu" "coo-sorting-gpu" "csr-gpu" "blcsr" "nwcsr")

error_exit() {
    echo -e "\e[31m$1\e[0m"
    exit 1
}

contains() {
    local element="$1"; shift
    for e in "$@"; do [[ "$e" == "$element" ]] && return 0; done
    return 1
}

prepare_output_file() {
    local res_file="$1"; shift
    local formats=("$@")
    echo -e "#nb-elt\t${formats[*]// /\\t}" > "$res_file"
}

replace_placeholders() {
    local template="$1" mesh="$2" output="$3"
    sed "s|MESH_FILE|${mesh}|g" "$template" > "$output"
}

append_formats() {
    local file="$1"; shift
    local formats=("$@")
    for format in "${formats[@]}"; do
        sed -i "/<!-- FORMATS -->/a <${format}>true</$format>" "$file"
    done
}

run_test() {
    local executable="$1" test_file="$2" mpi_num="$3" accelerated=$4 cache_warming="$5"
    local args="-A,CACHE_WARMING=${cache_warming}" # "-A,T=8" # for multiple threads
    $accelerated && args+=" -A,AcceleratorRuntime=${ACCELERATOR_RUNTIME}"
    if $IS_ADASTRA; then
      $accelerated && SRUN_ARGS+=" ${SRUN_ARGS_TO_ADD_IF_ACCELERATED}"
      srun --ntasks="$mpi_num" $SRUN_ARGS -- "$executable" "$test_file" $args > "stdout.txt" 2> "stderr.txt"
    else
      mpirun -n "$mpi_num" "$executable" "$test_file" $args > "stdout.txt" 2> "stderr.txt"
    fi
}

process_results() {
    local brief_file="$1" res_file="$2"; shift; shift
    formats=("$@")
    local line=$(grep "Element" "$brief_file" | awk '{print $2}')
    for format in "${formats[@]}"; do
        local time=$(grep "AssembleBilinearOperator_${format}:" "$brief_file" | awk '{print $2}')
        line+="\t${time:-NaN}"
    done
    echo -e "$line" >> "$res_file"
}

clear_output() {
    rm -rf "output/courbes" "output/depouillement" "fatal_4"
}

unique_directory_name() {
    local base_name="$1"
    local dir_name="$base_name"
    local counter=1

    # Check if the directory exists, and increment the counter until it doesn't
    while [[ -d "$dir_name" ]]; do
        dir_name="${base_name}.${counter}"
        ((counter++))
    done

    echo "$dir_name"
}

WORKING_DIR=$(unique_directory_name "$WORKING_DIR_BASE")

log_configuration() {
    local log_file="configuration.log"
    echo "Logging benchmark configuration to $log_file"
    {
        echo "============================================"
        echo "Benchmark Configuration"
        echo "============================================"
        echo "Executable: $EXECUTABLE"
        echo "Python Script: $PYTHON_SCRIPT"
        echo "Cache Warming: $CACHE_WARMING"
        echo "Working Directory: $WORKING_DIR"
        echo "Accelerator Runtime: $ACCELERATOR_RUNTIME"
        echo "MPI Configurations: ${MPI_N[*]}"
        echo "MPI Accelerated Configurations: ${MPI_N_ACCELERATED[*]}"
        echo "CPU Formats: ${CPU_FORMATS[*]}"
        echo "GPU Formats: ${GPU_FORMATS[*]}"
        echo "Dimensions: ${DIMENSIONS[*]}"
        echo "Mesh Sizes: ${SIZES[*]}"
        echo "Adastra Configuration: $IS_ADASTRA"
        echo "CPU Per Task: $CPU_PER_TASK"
        echo "Threads Per Core: $THREADS_PER_CORE"
        echo "\`srun\` Args: $SRUN_ARGS"
        echo "Additional \`srun\` Args if Accelerated: $SRUN_ARGS_TO_ADD_IF_ACCELERATED"
        echo "Templates:"
        for template in "${TEMPLATE_FILENAMES[@]}"; do
            echo "  - $template"
        done
        echo "2D Meshes:"
        for size in "${!MESHES_2D[@]}"; do
            echo "  $size: ${MESHES_2D[$size]}"
        done
        echo "3D Meshes:"
        for size in "${!MESHES_3D[@]}"; do
            echo "  $size: ${MESHES_3D[$size]}"
        done
        echo "============================================"
    } | tee "$log_file"
}

#======================================================================================
# Benchmarking Logic
#======================================================================================

mkdir -p "$WORKING_DIR" && cd "$WORKING_DIR" || error_exit "Failed to access working directory."

log_configuration

for dim in "${DIMENSIONS[@]}"; do
  echo "Processing ${dim}D meshes..."
  dim_dir="${dim}D"
  mkdir -p "$dim_dir" && cd "$dim_dir" || error_exit "Failed to enter $dim_dir."

  template_var="TEMPLATE_FILENAMES[$((dim-2))]"
  meshes_var="MESHES_${dim}D"

  for accelerated in false true; do

    if $accelerated; then
      mpi_array=("${MPI_N_ACCELERATED[@]}")
      format_array=("${GPU_FORMATS[@]}")
      all_format_array=("${ALL_GPU_FORMATS[@]}")
    else
      mpi_array=("${MPI_N[@]}")
      format_array=("${CPU_FORMATS[@]}" "${GPU_FORMATS[@]}")
      all_format_array=("${ALL_CPU_FORMATS[@]}" "${ALL_GPU_FORMATS[@]}")
    fi

    for mpi_n in "${mpi_array[@]}"; do
      dir="${mpi_n}-mpi-instance"
      $accelerated && dir="${ACCELERATOR_RUNTIME}.${dir}"
      mkdir -p "$dir" && cd "$dir" || error_exit "Failed to enter $dir."

      res_file="results.tsv"
      prepare_output_file "$res_file" "${all_format_array[@]}"

      for size in "${SIZES[@]}"; do
        mesh_file=$(eval "echo \${${meshes_var}[$size]}")

        [[ -e "$mesh_file" ]] || error_exit "Mesh file $mesh_file not found."

        $accelerated &&  test_name="Test.${dim}D.${ACCELERATOR_RUNTIME}.${mpi_n}-mpi-instance.${size}" || test_name="Test.${dim}D.${mpi_n}-mpi-instance.${size}"
        mkdir -p "$test_name" && cd "$test_name" || error_exit "Failed to enter $test_name."

        test_file="${test_name}.arc"
        replace_placeholders "${!template_var}" "$mesh_file" "$test_file"
        append_formats "$test_file" "${format_array[@]}"

        echo -n "Starting ${test_file}... "
        run_test "$EXECUTABLE" "$test_file" "$mpi_n" $accelerated "$CACHE_WARMING" || error_exit "\nTest $test_file failed."
        mv "./output/listing/time_stats.json" "./"
        python "$PYTHON_SCRIPT" "./time_stats.json" "BuildMatrix,AddAndCompute" > "brief.txt" || error_exit "An error occured in ${PYTHON_SCRIPT}."
        process_results "./brief.txt" "../$res_file" "${all_format_array[@]}"
        clear_output
        echo "Done"

        cd "../" || error_exit "Failed to return to size directory."
      done
      cd "../" || error_exit "Failed to return to mpi-instance directory."
    done
  done
  cd "../" || error_exit "Failed to return to dimension directory."
done

