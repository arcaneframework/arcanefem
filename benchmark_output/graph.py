import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
from matplotlib.colors import hsv_to_rgb

# size > method > samples

files = []
nb_simulations = 3
nb_sizes = 4
nb_attempts = 3

for i in range(nb_sizes):
	cur_size = []

	cur_size.append(f'petsc_direct_amgx_cuda_r{i}.out')
	cur_size.append(f'petsc_direct_gamg_cuda_r{i}.out')
	cur_size.append(f'petsc_direct_gamg_mpi_16_r{i}.out')
	cur_size.append(f'petsc_direct_gamg_mpi_20_r{i}.out')
	cur_size.append(f'hypre_direct_amg_mpi_16_r{i}.out')
	cur_size.append(f'hypre_direct_amg_mpi_20_r{i}.out')

	files.append(cur_size)

# one graph per size

# x-coordinate: method name
# y-coordinate: time spent

# colored bar graph with each bar having three colors:
# one for matrix
# one for vectors
# one for solve

sizes = []

for size in range(nb_sizes):
	cur_size = []

	for file in files[size]:
		# for each method, record matrix, vector and solve
		matrix = [[] for _ in range(nb_attempts)]
		vector = [[] for _ in range(nb_attempts)]
		solve = [[] for _ in range(nb_attempts)]

		for i in range(1, nb_simulations + 1):
			filepath = f'./benchmark_output{i}/{file}.csv'

			with open(filepath, newline='') as csvfile:
				reader = csv.DictReader(csvfile, delimiter=';')
				# print(f'in {filepath}: ')
				# only first attempt for now

				cur_attempt = 0

				for row in reader:
					matrix[cur_attempt].append(float(row['matrix_creation_time']))
					vector[cur_attempt].append(float(row['vector_creation_time']))
					solve[cur_attempt].append(float(row['solve_time']))
					cur_attempt += 1

		matrix = np.array(matrix)
		vector = np.array(vector)
		solve = np.array(solve)

		cur_method = np.array([matrix, vector, solve])
		cur_method = np.average(cur_method, axis=2)
		cur_size.append(cur_method)

		# cur_method[0] = [avg_mat_attempt1, avg_mat_attempt2, avg_mat_attempt3]
		# cur_method[1] = [avg_vec_attempt1, avg_vec_attempt2, avg_vec_attempt3]
		# cur_method[2] = [avg_solve_attempt1, avg_solve_attempt2, avg_solve_attempt3]

	sizes.append(cur_size)

cur_size_i = int(sys.argv[1])
sizes[cur_size_i] = np.nan_to_num(sizes[cur_size_i])
a = np.array(sizes[cur_size_i]).T
# a[attempt_i][element_i] = [method_1_element_i, method_2_element_i, ..., method_n_element_i]

part_names = ['Matrix', 'Vector', 'Solve']
attempt_names = [f'Attempt {i+1}' for i in range(a.shape[0])]
method_names = files[cur_size_i]

nb_parts = len(part_names)
nb_method = len(method_names)

bar_width = 0.25
x = np.arange(nb_method)

fig, ax = plt.subplots(figsize=(10, 6))
width = 0.2

# --- Create the 2D color gradient: green → red (H), darker → brighter (V)
H = np.linspace(0, 0.4, nb_parts)   # green (0.33) → red (0)
V = np.linspace(1.0, 0.5, nb_attempts) # darker → brighter per attempt

# Function to get color for part j, attempt i
def get_color(i, j):
    return hsv_to_rgb([H[j], 0.8, V[i]])

for attempt_i in range(nb_attempts):
    x_shifted = x + (attempt_i - nb_attempts / 2) * bar_width + bar_width / 2
    bottom = np.zeros(nb_method)
    for part_i in range(nb_parts):
        ax.bar(x_shifted, a[attempt_i, part_i, :],
               bar_width, bottom=bottom,
               color=get_color(attempt_i, part_i),
               label=f'{part_names[part_i]}' if attempt_i == 0 else "")
        bottom += a[attempt_i, part_i, :]

ax.set_xticks(x)
ax.set_xticklabels(method_names)
ax.set_ylabel("Execution time (s)")
ax.set_title("Time per method")
ax.legend(loc='upper left', fontsize='small')

plt.tight_layout()
plt.show()
