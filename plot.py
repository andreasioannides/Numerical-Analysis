import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from math import pi

def values(t: sp.Expr, x_t: sp.Expr, y_t: sp.Expr) -> tuple:
    t_values = np.linspace(0, 1, 1000)  # t ∈ [0, 1]
    x_values = np.zeros_like(t_values)
    y_values = np.zeros_like(t_values)

    for i, value in enumerate(t_values):
        x_values[i] = x_t.subs(t, value)
        y_values[i] = y_t.subs(t, value)

    return t_values, x_values, y_values

def plot_2d_bezier(cp: np.ndarray, t: sp.Expr, x_t: sp.Expr, y_t: sp.Expr):
    '''Plot the Bezier curve.'''

    t_values, x_values, y_values = values(t, x_t, y_t)

    plt.plot(x_values, y_values, color='red', linewidth=2.5)
    plt.scatter(cp[:, 0], cp[:, 1], color='black', s=50, facecolors='none', edgecolors='black')

    for i in range(cp.shape[0]-1):
        plt.plot([cp[i, 0], cp[i+1, 0]], [cp[i, 1], cp[i+1, 1]], "black", linewidth=1.5)
    
    plt.title('Bezier curve')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()

def plot_3d(t: sp.Expr, x_t: sp.Expr, y_t: sp.Expr):
    '''Plot the 3d revolution of the Bezier curve around y-axis.'''

    t_values, x_values, y_values = values(t, x_t, y_t)
    u = np.arange(0, 2*pi, 2*pi/360)

    xn = np.outer(x_values, np.cos(u))
    yn = np.outer(y_values, np.ones_like(u))
    zn = np.outer(x_values, np.sin(u))
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(xn, yn, zn, cmap='viridis')
    plt.plot(x_values, y_values, color='red', linewidth=3)

    ax.set_title('Revolved Bezier curve')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

def plot_leastSquares(nodes: np.ndarray, symbol: sp.Symbol, f: sp.Expr):  

    V_values = np.arange(1.5, 35, 0.5)
    m_values = np.zeros_like(V_values, dtype=float)

    for i, v in enumerate(V_values):
        result = f.subs(symbol, v).evalf()
        if sp.im(result) == 0:  # check if the result is real
            m_values[i] = np.real(result)
        else:
            m_values[i] = np.nan  # assign NaN for complex values

    plt.plot(V_values, m_values, color='red', linewidth=2)
    plt.scatter(nodes[:, 0], nodes[:, 1], color='black', s=20)
    
    plt.title('Least Squares')
    plt.xlabel('V')
    plt.ylabel('m')
    plt.grid(True)
    plt.show()

def plot_thomas(x: np.ndarray):
    x = np.flip(x)
    y = np.full(x.shape[0], 0)

    plt.vlines(0, 0, 70, linestyles='dashed', zorder=1)
    plt.scatter(y, x, color='red', edgecolor='black', s=90, zorder=2)

    labels = ['x12', 'x11', 'x10', 'x9', 'x8', 'x7', 'x6', 'x5', 'x4', 'x3', 'x2', 'x1']
    
    for i in range(int(len(labels)/2)):
        plt.text(y[2*i]+0.15, x[2*i], labels[2*i], va='center')

    for i in range(int(len(labels)/2)):
        plt.text(y[2*i+1]-0.30, x[2*i+1], labels[2*i+1], 
        va='center')
    
    plt.ylabel('x')
    plt.grid(True)
    plt.xticks(color='w')
    plt.ylim(x[0]+2, 0)
    plt.xlim(-2, 2)
    plt.show()

def plot_delta(sol: np.ndarray):
    v = [-1.1, 0.1, 1.1]
    m = 0 # max value of x

    for j, x in enumerate(sol):
        if (max(x) > m):
            m = max(x)

        x = np.flip(x)
        y = np.full(x.shape[0], v[j])

        plt.vlines(v[j], 0, x[0], linestyles='dashed', zorder=1)
        plt.vlines(v[j]+0.15, x[9], x[2], linestyles='dashed', zorder=1)
        plt.hlines(x[2], v[j], v[j]+0.15, linestyles='dashed', zorder=1)
        plt.hlines(x[9], v[j], v[j]+0.15, linestyles='dashed', zorder=1)
        plt.scatter(y, x, color='red', edgecolor='black', s=90, zorder=2)

        labels = ['x12', 'x11', 'x10', 'x9', 'x8', 'x7', 'x6', 'x5', 'x4', 'x3', 'x2', 'x1']
        
        for i in range(int(len(labels)/2)):
            if (i == 1):
                plt.text(y[2*i]-0.4, x[2*i], labels[2*i], va='center') # plot x10

        for i in range(int(len(labels)/2)):
            if (i == 4):
                plt.text(y[2*i+1]-0.3, x[2*i+1], labels[2*i+1], va='center') # plot x3

    plt.ylabel('x')
    plt.grid(True)
    plt.xticks(color='w')
    plt.ylim(m+2, 0)
    plt.xlim(-2, 2)
    plt.show()