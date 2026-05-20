import numpy as np
import matplotlib.pyplot as plt

# Generate x values
x = np.linspace(-5, 5, 500)

# Gaussian parameters
mu = 0       # mean
sigma = 1    # standard deviation
A = 1        # amplitude

# Gaussian function
y = A * np.exp(-((x - mu)**2) / (2 * sigma**2))

# Plot
plt.plot(x, y, linewidth=2)
plt.title("Gaussian Curve")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
plt.show()