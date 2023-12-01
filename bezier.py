import numpy as np
from math import factorial as fact
import sympy as sp
from math import pi
from numpy.polynomial.legendre import leggauss

def factorial(n: int) -> int:  # instead use math.factorial
    '''Calculate the factorial of a given number.'''

    if (n < 0):
        raise ValueError("Factorial is not defined for negative numbers.")
    elif (n == 0):
        return 1
    else:
        fact = n
        for i in range(1, n):
            fact = fact*(n-i)
    
        return fact

def binomial_coeff(n: int, i: int) -> int:
    '''Calculate the binomial coefficient given two numbers.'''

    n_fact = factorial(n)
    i_fact = factorial(i)    
    ni_fact = factorial(n-i)

    return n_fact/(i_fact*ni_fact)

def m_matrix(N: int) -> np.ndarray:
    '''Find the (N+1)x(N+1) matrix with coefficients which are used to find the C matrix.'''

    mtrx = np.zeros(shape=[N+1, N+1], dtype=np.int16)

    for i in range(N+1):
        for j in range(N+1):
            if (j >= i):
                mtrx[i, j] = pow(-1, j-i) * binomial_coeff(N,j) * binomial_coeff(j,i)

    return mtrx

def t_vector(N: int) -> np.ndarray:
    '''Initialize the vector with the parameter t.'''

    t_vector = np.empty(shape=[N+1, 1], dtype=object) # vector with the parameter t raised to i-power

    t = sp.symbols('t')

    for i in range(N+1):
        t_vector[i] = t**i
    
    return t, t_vector

def parametric_equations(m: np.ndarray, t_vector: np.ndarray, cp: np.ndarray) -> sp.Expr:
    '''Find the parametric equation x(t) or y(t) of a Bezier curve given the coordinate x or y of N+1 control points (CP).'''

    C_matrix = np.matmul(m, t_vector)

    param = 0 # parametric equation

    for i in range(C_matrix.shape[0]):
        param = param + C_matrix[i] * cp[i]

    return param[0]  # return the param as a single expression. Note that the param is an array

def bezier(N: int, cp: np.ndarray) -> tuple:
    mtrx = m_matrix(N)
    t, t_vec = t_vector(N)

    x_t = parametric_equations(mtrx, t_vec, cp[:, 0])
    y_t = parametric_equations(mtrx, t_vec, cp[:, 1])

    print("Parametric equations:")
    print(f"x(t)={x_t}")
    print(f"y(t)={y_t}\n")

    return t, x_t, y_t

def simpson_38(a: float, b: float, t: sp.Expr, x_t: sp.Expr, y_t: sp.Expr, n: int) -> float:
    '''Find the integral of a function f within [a,b] with Simpson 3/8 method.'''

    f = pi * abs(pow(x_t, 2) * sp.diff(y_t, t))
    h = (b - a) / n
    
    t_values = np.linspace(a, b, n)
    f_values = np.zeros_like(t_values)

    for i, value in enumerate(t_values):
        f_values[i] = f.subs(t, value)

    s1 = sum(f_values[1:n:3])
    s2 = sum(f_values[2:n:3])
    s3 = sum(f_values[3:n-1:3])

    I = (3*h)/8 * (f.subs(t, t_values[0]) + 3*s1 + 3*s2 + 2*s3 + f.subs(t, t_values[n-1]))

    return I

def gauss_integral(x_t: sp.Expr, y_t: sp.Expr, t: sp.Expr, n: int) -> np.ndarray:
    '''Find the integral of a function f with the Gauss integration method.'''

    f = pi * abs(pow(x_t, 2) * sp.diff(y_t, t))
    td = sp.symbols('td')
    f = f.subs(t, 0.5+0.5*td) * 0.5
    x, w = leggauss(n) # calculate the integration points and associated weights
    I = 0

    for x_i, w_i in zip(x, w):
        I += w_i * f.subs(td, x_i)

    return I