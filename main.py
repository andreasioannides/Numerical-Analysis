from bezier import bezier, simpson_38, gauss_integral
from plot import plot_2d_bezier, plot_3d, plot_leastSquares, plot_thomas, plot_delta
from leastSquare import least_squares, val
from springs import thomas, delta, tridiagonal
import numpy as np
import sympy as sp

def main():
    N = 6
    cp = np.array([[0, 10], [1.8, 8.01], [2.9, 7], [8, 4.7], [4.3, 4.7], [1.9, 1.9], [0, 0.5]])

    # There are 6 problems. Unhighlight the one you want to solve.

    '''Problem A: Bezier curve'''
    t, x_t, y_t = bezier(N, cp)

    plot_2d_bezier(cp, t, x_t, y_t)
    plot_3d(t, x_t, y_t)

    a = 0
    b = 1
    nodes = [7, 55, 208, 508]

    '''Problem B: Volume'''
    # print("Simpson 3/8:")
    # for n in nodes:
    #     I = simpson_38(a, b, t, x_t, y_t, n)
    #     print(f"{n} nodes: I={I}.")

    '''Problem C: Least Squares'''
    # nodes = np.array([[2, 2.1], [4.8, 2.16], [7.6, 2.39], [10.4, 3.25], [13.2, 4.65], [16, 6.59], [18.8, 9.3], [21.6, 12.44], [24.4, 16.6], [27.2, 21.29], [30, 26.7]])
    # iters = 400
    # epsilon = 1.0e-10

    # c1, c2, c3, c4, c5 = sp.symbols('c1 c2 c3 c4 c5')
    # V = sp.symbols('V')
    # m = sp.symbols('m')

    # symbols = np.array([c1, c2, c3, c4, c5, V, m], dtype=object)
    # g = c1 + (c2*pow(V, 2) + c3*V + c4)*sp.log(c5*V, 10)

    # # initial = np.array([0.1, -0.1, 0.2, -0.2, 0.3]) 
    # initial = np.array([1.1, 0.9, -1.1, -0.9, 1])
    # # initial = np.array([10, 5, -3, 0, 0.05])

    # c_values = least_squares(nodes, symbols, g, initial, iters, epsilon)

    # # np.savetxt('c1.txt', c_values) # save c to a txt file
    # # c_values = np.loadtxt('c1.txt') # load the txt file

    # values = dict(zip(symbols[:-2], c_values))
    # g_v = val(g, values)
    # print(c_values)
    # print(val(g_v, {V: 26.45624}))
    
    # plot_leastSquares(nodes, symbols[5], g_v)

    ''''Problem D: Thomas algorithm'''
    # n = 12
    # k = 218.3326

    # b = np.zeros(shape=(n)) # main diagonal
    # for j in range(n):
    #     b[j] = 2*k

    # b[-1] = k

    # a = np.zeros(shape=(n-1)) # lower diagonal
    # for j in range(n-1):
    #     a[j] = -k

    # mass = 20.0305
    # grav = 9.81
    
    # d = np.zeros(shape=(n))
    # for j in range(n):
    #     if (j == n-1):
    #         d[n-1] = mass*grav + k
    #     else:
    #         d[j] = mass * grav

    # x = thomas(a, b, a, d)
    # print(x)
    # plot_thomas(x)

    '''Problem E: Delta formulation'''
    # n = 12
    # k = 218.3326
    # k_star_list = [50, k, 300]
    # initial = np.array([11.8, 22.7, 32.7, 41.8, 50, 57.3, 63.7, 69.2, 73.8, 77.5, 80.3, 82.2])
    # iters = 3000
    # epsilon = 1.0e-10
    # mass = 20.0305
    # grav = 9.81
    # sols = np.zeros((3, 12)) # [x1,...,x12] for each k*

    # for i, k_star in enumerate(k_star_list):
    #     # main diagonal
    #     b = np.zeros(shape=(n)) 
    #     for j in range(n):
    #         if (j == 2 or j == 9):
    #             b[j] = 2*k + k_star
    #         elif (j == n-1):
    #             b[j] = k
    #         else:
    #             b[j] = 2*k

    #     # lower diagonal
    #     a = np.zeros(shape=(n-1)) 
    #     for j in range(n-1):
    #         a[j] = -k
        
    #     d = np.zeros(shape=(n))
    #     for j in range(n):
    #         if (j == 2):
    #             d[2] = mass*grav - 7*k_star 
    #         elif (j == 9):
    #             d[9] = mass*grav + 7*k_star 
    #         elif (j == 11):
    #             d[11] = mass*grav + k
    #         else:
    #             d[j] = mass*grav

    #     A = tridiagonal(n, b, a, a)
    #     A[2, 9] = -k_star
    #     A[9, 2] = -k_star

    #     x = delta(A, d, initial, iters, epsilon)
    #     print(x)
    #     print()
    #     sols[i] = x

    # plot_delta(sols)

    '''Problem F: Gauss Integral'''
    # for i in range(2, 8):
    #     I = gauss_integral(x_t, y_t, t, i)
    #     print(f"n={i}: I={I}")

if __name__ == '__main__':
    main()