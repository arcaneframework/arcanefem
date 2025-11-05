import matplotlib.pyplot as plt
import numpy as np
import csv
import sys

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

print(sys.argv)
cur_size_i = int(sys.argv[1])
element_idx = int(sys.argv[2])
elements = ('Matrix', 'Vector', 'Solve')

method_names = tuple(files[cur_size_i])

for i in range(len(sizes[cur_size_i])):
	size = sizes[cur_size_i][i]
	sizes[cur_size_i][i] = size[element_idx]

sizes[cur_size_i] = np.nan_to_num(sizes[cur_size_i])

a = np.array(sizes[cur_size_i]).T

time_means = dict()

for i in range(nb_attempts):
	time_means[f'Attempt {i + 1}'] = tuple(a[i])

x = np.arange(len(method_names))
width = 0.25
multiplier = 0

fig, ax = plt.subplots(layout='constrained')

for attribute, measurement in time_means.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute)
    ax.bar_label(rects, padding=3)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Time (s)')
ax.set_title(f'{elements[element_idx]} time per method')
ax.set_xticks(x + width, method_names)
ax.legend(loc='upper left')
ax.set_ylim(0, np.max(a) * 1.05)

plt.show()

