#!/bin/bash

error_exit() { echo -e "\e[31m$1\e[0m"; exit 1; }
usage() { echo "Usage: $0 <file>"; exit 1; }

[ "$#" -eq 1 ] || usage
file=$1
[[ -f "$file" ]] || error_exit "Error: '${file}' is not a valid file"

generate_plot_title() {
    local filepath="$1"
    benchmark=$(echo "$filepath" | awk -F'/' '{print $(NF-3)}')
    filename=$(basename "$filepath")
    num_mpi_instances=$(echo "$filepath" | grep -oP '\d+(?=-mpi-instance)')
    num_cache_warming=$(echo "$filepath" | grep -oP '\d+(?=-cw)')
    dimension=$(echo "$filepath" | grep -oP '(?<=/)(2D|3D)(?=/)')
    runtime="cpu" # Default runtime
    [[ "$filepath" =~ cuda|hip ]] && runtime="gpu"
    [[ "$runtime" =~ gpu && "$filepath" =~ cuda ]] && runtime+=" (cuda)"
    [[ "$runtime" =~ gpu && "$filepath" =~ hip ]] && runtime+=" (hip)"
    title="$benchmark | $dimension, ${num_mpi_instances} MPI instance(s), ${num_cache_warming} cache warming(s), Runtime: $runtime"
}

is_col_empty() {
    col="$1"
    nan_count=$(awk -v col="$col" '{ if ($col == "NaN") count++ } END { print count }' "$file")
    [[ "$nan_count" -eq "$((num_rows - 1))" ]] && return 0 || return 1
}

generate_plot() {
    local mode="$1" col_end="$2" title_arr=("${!3}") pt_arr=("${!4}") dt_arr=("${!5}")
    echo "Info: Generating $mode plot..."
    cmd="set title '${title}' \nset xlabel '# Elements' \nset ylabel 'Execution Time (s)' \nset grid\n"
    cmd+="set logscale x \nset logscale y \nset key outside right \nplot "
    for i in $(seq 2 "$col_end"); do
        if ! is_col_empty "$i"; then
            j=$((i - 2))
            cmd+="'${file}' using 1:${i} with linespoints lw 2 pt ${pt_arr[$j]} dt ${dt_arr[$j]} title '${title_arr[$j]}',"
        fi
    done
    echo -e "${cmd%,}" | gnuplot -persist
    echo "Info: Done generating $mode plot"
}

generate_plot_title "$file"
cpu_flag=true gpu_flag=false
[[ "$runtime" =~ gpu ]] && cpu_flag=false gpu_flag=true

num_columns=$(head -n 1 "$file" | awk '{ print NF }')
num_rows=$(awk 'END { print NR }' "$file")
echo -e "Info: ${num_columns} columns in '$file'"
echo -e "Info: ${num_rows} rows in '$file'"

if $cpu_flag; then
    cpu_titles=("DOK" "COO_{CT}^{C}" "S-COO_{CT}^{C}" "CSR_{CT}^{C}" "COO_{CTG}^{C}" "S-COO_{CTG}^{C}" "CSR_{CTG}^{C}" "BL-CSR_{CTG}^{C}" "NW-CSR_{CTG}^{C}")
    cpu_pts=(1 2 3 4 2 3 4 5 6)
    cpu_dts=(1 1 1 1 2 2 2 2 2)
    generate_plot "CPU" 10 cpu_titles[@] cpu_pts[@] cpu_dts[@]
fi

if $gpu_flag; then
    gpu_titles=("COO_{CTG}^{G}" "S-COO_{CTG}^{G}" "CSR_{CTG}^{G}" "BL-CSR_{CTG}^{G}" "NW-CSR_{CTG}^{G}")
    gpu_pts=(2 3 4 5 6)
    gpu_dts=(1 1 1 1 1)
    generate_plot "GPU" 6 gpu_titles[@] gpu_pts[@] gpu_dts[@]
fi
