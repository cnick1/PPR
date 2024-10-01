import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches
from scipy.integrate import odeint
from examples.example11ExportedControlFunction import u
from sympy import *

g = 9.81
l = 10
m = 1
theta0 = -np.pi
xmax = 3.5
x1_range = np.linspace(-xmax, xmax, 400)
x2_range = np.linspace(-5, 5, 200)
X1, X2 = np.meshgrid(x1_range, x2_range)
U, V = np.zeros_like(X1), np.zeros_like(X2)

def f(x1, x2):
    return [x2, 3 * u(x1, x2) / m / l**2 + 3 * g / 2 / l * np.sin(x1)]

def pendulum_derivatives(x, t):
    x1, x2 = x
    x1_dot, x2_dot = f(x1, x2)
    return [x1_dot, x2_dot]

def plot_function(x1_range, x2_range, f, xmax, theta0, g, l):
    U = np.zeros((len(x2_range), len(x1_range)))
    V = np.zeros((len(x2_range), len(x1_range)))
    for i in range(len(x1_range)):
        for j in range(len(x2_range)):
            x1 = x1_range[i]
            x2 = x2_range[j]
            U[j, i], V[j, i] = f(x1, x2)

    fig, ax = plt.subplots()
    ax.streamplot(X1, X2, U, V, color='lightgray', density=.8, linewidth=1, cmap=plt.cm.gray,broken_streamlines=False)
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(-5, 5)

    x0 = [theta0, np.sqrt(3 * g * (1 - np.cos(theta0)) / l)]
    x_t = odeint(pendulum_derivatives, x0, np.linspace(0, 5, 1000))
    ax.plot(x_t[:, 0], x_t[:, 1], 'r-', linewidth=2)
    
    ax.set_frame_on(False)
    plt.axis('off') # this rows the rectangular frame 
    ax.get_xaxis().set_visible(False) # this removes the ticks and numbers for x axis
    ax.get_yaxis().set_visible(False) # this removes the ticks and numbers for y axis

    for art in ax.get_children():
        if not isinstance(art, matplotlib.patches.FancyArrowPatch):
            continue
        art.remove()        # Method 1

plot_function(x1_range, x2_range, f, xmax, theta0, g, l)
plt.savefig(f"plots/example11_phasePortrait_u{int(controlDegree)}.pdf", bbox_inches='tight', pad_inches=0, transparent=True)
plt.savefig("plots/example11_phasePortrait.pdf", bbox_inches='tight', pad_inches=0, transparent=True)
plt.close('all')
x1, x2 = symbols('x1 x2')
print("u(x1,x2) = " + str(u(x1,x2)))