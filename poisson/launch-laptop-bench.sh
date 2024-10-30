#!/bin/bash

working_dir="Benchmark"

# Must be in absolute patawk '/Loop/{getline; if(/Fem/) {print; exit}}' "$log_file"h
executable="$HOME/Git/toutane/arcanefem_forked/build/poisson/Poisson"
template_filename="$HOME/Git/toutane/arcanefem_forked/poisson/TEST_TEMPLATE.xml"

mesh_2D_small="$HOME/Git/toutane/arcanefem_forked/meshes/msh/L-shape.msh"
mesh_2D_medium=""
mesh_2D_large=""

mesh_3D_small=""
mesh_3D_medium=""
mesh_3D_large=""

formats=("legacy" "coo" "coo-sorting" "csr" "blcsr" "nwcsr")
# formats=("coo")
# sizes=("small" "medium" "large")
sizes=("small")
# cpu_core_numbers=(1 2 5 10 20)
cpu_core_numbers=(1)

[[ ! -d "$working_dir" ]] && mkdir "$working_dir"
cd "$working_dir"

launchCpuTests () # arg1: test filename
{
  for cpu_n in "${cpu_core_numbers[@]}"; do

    dir="cpu_$cpu_n"
    [[ ! -d "$dir" ]] && mkdir "$dir"
    cd "$dir"

    mpirun -n "$cpu_n" "$executable" "../$1" > stdout.txt 2> stderr.txt
    return_code=$(echo $?)

    if [ "$return_code" -eq 0 ]; then
      echo -e "$test_filename [$cpu_n CPU]: SUCCESS"

      log_file="./output/listing/logs.0"

      grep -E "^ *Fem" "$log_file" | tail -n 2 | head -n 1 |  awk '{print $1, $2}'
      grep -E "^ *Assemble.*BilinearOperator" "$log_file" | awk '{print $1, $2}'
      grep -E "^[ a-zA-Z]*BuildMatrix" "$log_file" | awk '{print $1, $2}'

      echo ""

    else
      echo -e "\e[31m$test_filename [$cpu_n CPU]: FAILURE"
      echo -e "command: mpirun -n $cpu_n $executable ../$1\e[0m"
    fi

    cd "../"

  done
}


# for dim in {2..3}; do
for dim in 2; do

  dir="${dim}D"
  [[ ! -d "$dir" ]] && mkdir "$dir"
  cd "$dir"

  for size in "${sizes[@]}"; do

    dir="$size"
    [[ ! -d "$dir" ]] && mkdir "$dir"
    cd "$dir"

    for format in "${formats[@]}"; do

      dir="$format"
      [[ ! -d "$dir" ]] && mkdir "$dir"
      cd "$dir"

      test_name="Test.${dim}D.${size}.${format}"
      test_filename="${test_name}.arc"

      mesh_file="mesh_${dim}D_${size}"

      sed -e "s|MESH_FILE|${!mesh_file}|g" -e "s|FORMAT|${format}|g" "$template_filename" > "$test_filename";

      launchCpuTests $test_filename

      cd "../"
    done
    cd "../"
  done
  cd "../"
done
