import numpy as np
import matplotlib.pyplot as plt

# Example data (same as before)
a = np.array([
    [[2.22121239e-01, 2.23719438e-01, 1.12501621e-01, 1.36573553e-01, 5.84764481e-02, 7.04135100e-02],
     [3.84017626e-02, 3.76919905e-02, 3.80110741e-03, 6.66149457e-03, 1.70381864e-03, 5.69407145e-03],
     [4.48013624e-01, 1.39632686e+00, 7.82595793e-01, 1.05990052e+00, 1.47913059e-01, 1.62120263e-01]],

    [[1.88413223e-01, 1.88172976e-01, 1.08112176e-01, 1.27811750e-01, 5.75365225e-02, 6.96541468e-02],
     [3.69394620e-02, 3.73817285e-02, 3.87899081e-03, 4.35090065e-03, 8.35895538e-04, 3.23994954e-03],
     [3.82056952e-01, 1.33728091e+00, 7.79146910e-01, 9.29013968e-01, 1.27653519e-01, 1.42678420e-01]],

    [[1.87301636e-01, 1.79371039e-01, 1.09248718e-01, 1.38278087e-01, 5.82503478e-02, 5.51739534e-02],
     [3.68132591e-02, 3.72581482e-02, 3.88153394e-03, 4.20729319e-03, 8.31683477e-04, 6.57359759e-03],
     [3.82093191e-01, 1.33739193e+00, 7.69998550e-01, 8.50647449e-01, 1.27562920e-01, 1.26240810e-01]]
])

# Labels (example)
part_names = ['Matrix', 'Vector', 'Solve']
attempt_names = [f'Attempt {i+1}' for i in range(a.shape[0])]
method_names = [f'Method {i+1}' for i in range(a.shape[2])]

# Dimensions
n_attempts, n_parts, n_methods = a.shape
bar_width = 0.25
x = np.arange(n_methods)

# Color scheme (one color per method)
colors = plt.cm.tab10(np.linspace(0, 1, n_parts))

fig, ax = plt.subplots(figsize=(10, 6))

# Plot
for attempt_i in range(n_attempts):
    x_shifted = x + (attempt_i - n_attempts / 2) * bar_width + bar_width / 2
    bottom = np.zeros(n_methods)
    for method_i in range(n_parts):
        ax.bar(x_shifted, a[attempt_i, method_i, :],
               bar_width, bottom=bottom,
               color='#ff0000',
               label=f'{part_names[method_i]}' if attempt_i == 0 else "")
        bottom += a[attempt_i, method_i, :]

# Styling
ax.set_xticks(x)
ax.set_xticklabels(method_names)
ax.set_ylabel("Execution time (s)")
ax.set_xlabel("Program parts")
ax.set_title("Benchmark results by program part, attempt, and method")
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')

plt.tight_layout()
plt.show()
