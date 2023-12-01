import numpy as np
import sympy as sp

def val(f: sp.Expr, value: dict) -> any:
    '''Assign a value to an expression.'''

    return f.subs(value)

def der(f: sp.Expr, symbol: sp.Expr) -> any:
    '''Calculate the derivative of a function f respect to its variable.'''

    return sp.diff(f, symbol)

def Crout(A: np.ndarray) -> tuple:
    '''LU decomposition with Crout method (decomposition of A to L and U matrices).'''

    n = A.shape[0]

    L = np.zeros_like(A, dtype=np.float64)
    U = np.zeros_like(A, dtype=np.float64)
    
    for i in range(n):
        U[i, i] = 1

    for i in range(n):
        L[i,0] = A[i,0]

    for j in range(1, n):
        U[0,j] = A[0,j] / L[0,0]

    for r in range(1, n):
        for i in range(r, n):
            s1 = sum(L[i,k] * U[k,r] for k in range(0, r))
            L[i,r] = A[i,r] - s1

        for j in range(r+1, n):
            s2 = sum(L[r,k] * U[k,j] for k in range(0, r))
            U[r,j] = (A[r,j] - s2) / L[r,r]

    return L, U

def forward_sub(L: np.ndarray, b: np.ndarray) -> np.ndarray:
    '''Forward substitution L*y=b in LU method.'''

    y = np.zeros(shape=b.shape) 

    y[0] = b[0] / L[0,0]

    for i in range(0, b.shape[0]):
            s = sum(L[i, j] * y[j] for j in range(0, i))
            y[i] = (1/L[i,i]) * (b[i] - s) 

    return y

def backward_sub(U: np.ndarray, y: np.ndarray) -> np.ndarray:
    '''Backward substitution U*x=y in LU method.'''

    n = y.shape[0]
    x = np.zeros(shape=n) 

    x[n-1] = y[n-1] / U[n-1,n-1]

    for i in range(n-2, -1, -1):
            s = sum(U[i, j] * x[j] for j in range(i+1, n))
            x[i] = (1/U[i,i]) * (y[i] - s) 

    return x

def LU(A: np.ndarray, b: np.ndarray, crout: bool =True) -> np.ndarray:
    '''LU factorization method solves systems of linear equations.
       By default Crout method is used to decompose A matrix.'''

    if (crout):
        L, U = Crout(A)
    # else: L, U =  doolittle(J, F)

    y = forward_sub(L, b)
    x = backward_sub(U, y)

    return x

def jacobian(symbols: np.ndarray, f: np.ndarray, c_old: np.ndarray) -> np.ndarray:
    '''Create the Jacobian matrix.'''

    J = np.zeros(shape=[f.shape[0], symbols.shape[0]])
    values = dict(zip(symbols, c_old))
    
    for i, f_i in enumerate(f):
        for j, c_j in enumerate(symbols):
            fp = der(f_i, c_j)  # partial derivative: Ji=∂fi/∂ci
            J[i, j] = val(fp, values)
    
    return J

def F_vector(symbols: np.ndarray, f: np.ndarray, c_old: np.ndarray) -> np.ndarray:
    '''Create the F vector with f[i](c1,...,cn) used to calculate the dc vector.'''

    F = np.zeros(shape=(f.shape[0]))
    substitutions = {symbol: c_value for symbol, c_value in zip(symbols, c_old)}

    for j in range(f.shape[0]):
        F[j] = val(f[j], substitutions)

    return F

def newraph(initial: np.ndarray, iters: int, epsilon: float, symbols: np.ndarray, f: np.ndarray) -> np.ndarray:
    '''Newton-Raphson method solves systems of nonlinear equations.'''

    c_old = initial
    conv = False

    for i in range(iters):
        print(f"Iter:{i}, c={c_old}")

        J = jacobian(symbols, f, c_old)
        F = F_vector(symbols, f, c_old)
        dc = LU(J, -F)
        c_new = c_old + dc

        if (abs((c_new - c_old)/c_new) <= epsilon).all(): # converged=True
            conv = True
            break
        
        c_old = c_new

    print(f"Converged: {conv}.")

    return c_new

def least_squares(nodes: np.ndarray, symbols: np.ndarray, g: sp.Expr, initial: np.ndarray, iters: int, epsilon: float) -> np.ndarray:
    '''Determine the best fitting a curve on top of given nodes with the Least Squares method.'''

    m = symbols[-1] # -1: access the last element of a numpy array which is m
    V = symbols[-2]
    
    e = g - m
    E = 0

    for i in range(nodes.shape[0]):
        E += val(e, {V: nodes[i, 0], m: nodes[i, 1]}) ** 2
    
    f = np.zeros(shape=(symbols.shape[0]-2), dtype=object) # partial derivatives of E with respect to its variables, ex: f1=∂E/∂c1

    for i in range(symbols.shape[0]-2):
        f[i] = der(E, symbols[i])
    
    c_values = newraph(initial, iters, epsilon, symbols[:-2], f)

    return c_values