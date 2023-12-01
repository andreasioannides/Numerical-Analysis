import numpy as np

def thomas(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> np.ndarray:
    '''Solve a tridiagonal system of equations using the Thomas algorithm. Parameters a, b and c refer to lower, main and upper diagonals respectively.'''
    
    n = d.shape[0]

    a = np.copy(a)
    b = np.copy(b)
    c = np.copy(c)
    d = np.copy(d)

    # Forward substitution
    c[0] = c[0] / b[0]  # Divide first line by b[0]
    d[0] = d[0] / b[0]
    b[0] = 1
    
    for i in range(1, n):
        try:
            r = a[i-1] / b[i-1]
            b[i] = b[i] - r * c[i-1]
            d[i] = d[i] - r * d[i-1]
        except ZeroDivisionError:
            print("Divide by zero.")
            return
    
    # Backward substitution
    d[n-1] = d[n-1] / b[n-1]
    
    for i in range(n-2, -1, -1):
        d[i] = (d[i] - c[i] * d[i+1]) / b[i]
        
    return d

def tridiagonal(n: int, main: np.ndarray, upper: np.ndarray, lower: np.ndarray) -> np.ndarray:
    '''Create a tridiagonal matrix.'''

    main_diag = np.full(n, main)
    upper_diag = np.full(n-1, upper)
    lower_diag = np.full(n-1, lower)
    tridiagonal_matrix = np.diag(main_diag) + np.diag(upper_diag, k=1) + np.diag(lower_diag, k=-1)
    return tridiagonal_matrix

def delta(A: np.ndarray, b: np.ndarray, initial: np.ndarray, iters: int, epsilon: float) -> np.ndarray:
    '''Delta formulation method.'''

    main = np.diag(A)
    upper = np.diag(A, k=1)
    lower = np.diag(A, k=-1)
    A_dot = tridiagonal(A.shape[0], main, upper, lower) # approximate matrix 

    x_old = initial
    conv = False

    for i in range(iters):
        # print(f"Iter:{i}, x={x_old}")

        d = b - np.dot(A, x_old)
        dx = thomas(lower, main, upper, d)
        x_new = x_old + dx

        if (abs((x_new - x_old)/x_new) <= epsilon).all():
            conv = True
            break
        
        x_old = x_new

    print(i)
    print(f"Converged: {conv}.")

    return x_new