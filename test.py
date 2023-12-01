import numpy as np
from numpy.polynomial.legendre import leggauss

def gauss_legendre_integration(n):
    # Calculate integration points and weights using leggauss
    x, w = leggauss(n)
    return x, w

# Example usage
n = 6 # Number of integration points (can be adjusted)
integration_points, weights = gauss_legendre_integration(n)

# Print the integration points and weights
print("Integration Points:", integration_points)
print("Weights:", weights)
