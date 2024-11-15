#!/bin/bash

usage() {
    echo "Usage: plot-results.sh <file> <mpi-n> <cache-warming> <dimension> [-cpu | -gpu]"
    exit 1
}

if [[ $# -lt 3 ]]; then
    echo -e "Error: Missing arguments"
    usage
fi

file=$1
if [[ ! -f "$file" ]]; then
  echo -e "Error: '${file}' is not a valid file"
  exit 1
fi
shift

mpiNum=$1
shift

cacheWarming=$1
shift

dimension=$1
shift

cpu_flag=false
gpu_flag=false
# Process the flag
while [[ $# -gt 0 ]]; do
    case "$1" in
        -cpu)
            cpu_flag=true
            ;;
        -gpu)
            gpu_flag=true
            ;;
        *)
            echo -e "Error: Invalid option '$1'"
            usage
            ;;
    esac
    shift  # Move to the next argument
done

if [[ -z cpu_flag || -z gpu_flag ]]; then
    echo -e "Error:\tYou must specify either -cpu or -gpu"
    usage
fi
echo -e "Info: File: $file"

num_columns=$(head -n 1 "$file" | awk '{ print NF }')
num_rows=$(awk -v col=1 'END { print NR }' "$file")

echo -e "Info: ${num_columns} columns in '$file'"
echo -e "Info: ${num_rows} rows in '$file'"

function isColEmpty() {
  col="$1"
  nan_count=$(awk -v col="$col" '{ if ($col == "NaN") count++ } END { print count }' "$file")

  if [[ "$nan_count" -eq "((${num_rows} - 1))" ]]; then
    echo -e "Info: Column $col (${title_arr[(($col - 2))]}) is empty (full of 'NaN'): it will not be plotted"
    return 0
  fi

  return 1
}

title_arr=("DOK" "COO_{CT}^{C}" "S-COO_{CT}^{C}" "CSR_{CT}^{C}" "COO_{CTG}^{C}" "S-COO_{CTG}^{C}" "CSR_{CTG}^{C}" "BL-CSR_{CTG}^{C}" "NW-CSR_{CTG}^{C}" "COO_{CTG}^{G}" "S-COO_{CTG}^{G}" "CSR_{CTG}^{G}" "BL-CSR_{CTG}^{G}" "NW-CSR_{CTG}^{G}")
pt_arr=(1 2 3 4 2 3 4 5 6 2 3 4 5 6)
dt_arr=(1 1 1 1 2 2 2 2 2 1 1 1 1 1)

if $cpu_flag; then
  echo "Info: Generating CPU plot..."

  cmd="set title 'Format Comparison (Dimension: ${dimension}, CPU runtime, ${mpiNum} MPI instance, ${cacheWarming} cache warming)' \n"
  cmd+="set terminal wxt size 800,600 \n"
  cmd+="set xlabel '# Elements' \n"
  cmd+="set ylabel 'Execution Time (s)' \n"
  cmd+="set logscale x \n"
  cmd+="set logscale y \n"
  cmd+="set key outside right \n"
  cmd+="plot "

  for i in $(seq 2 10); do 
    if ! isColEmpty "$i"; then
      j=$(($i - 2))
      cmd+="'${file}' using 1:${i} with linespoints lw 2 pt ${pt_arr[$j]} dt ${dt_arr[$j]} title '${title_arr[$j]}',"
    fi
  done

  echo -e "$cmd" | gnuplot -persit
  echo "Info: Done generating CPU plot"
fi

if $gpu_flag; then
  echo "Info: Generating GPU plot..."

  cmd="set title 'Format Comparison (Dimension: ${dimension}, GPU runtime, ${mpiNum} MPI instance, ${cacheWarming} cache warming)' \n"
  cmd+="set terminal wxt size 800,600 \n"
  cmd+="set xlabel '# Elements' \n"
  cmd+="set ylabel 'Execution Time (s)' \n"
  cmd+="set logscale x \n"
  cmd+="set logscale y \n"
  cmd+="set key outside right \n"
  cmd+="plot "

  for i in $(seq 11 15); do 
    if ! isColEmpty "$i"; then
      j=$(($i - 2))
      cmd+="'${file}' using 1:${i} with linespoints lw 2 pt ${pt_arr[$j]} dt ${dt_arr[$j]} title '${title_arr[$j]}',"
    fi
  done

  echo -e "$cmd" | gnuplot -persit
  echo "Info: Done generating GPU plot"
fi
