import matplotlib.pyplot as plt
import numpy as np
import csv

# size > method > samples

files = []

for i in range(4):
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

for size in range(4):
	cur_size = []
	for file in files[size]:
		# for each method, record matrix, vector and solve
		cur_method = []
		for i in range(1, 4):
			filepath = f'./benchmark_output{i}/{file}.csv'
			with open(filepath, newline='') as csvfile:
				reader = csv.reader(csvfile, delimiter=';')
				print(f'in {filepath}: ')
				attempt = 0
				for row in reader:
					row['matrix_creation_time']
					attempt += 1
		cur_size.append(cur_method)
	sizes.append(cur_size)
