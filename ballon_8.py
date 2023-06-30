from math import cos, sin, pi, degrees
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
from celluloid import Camera


# Some const and vars
TIMECONST = 2.5
STEP1 = 0.01
STEP2 = 0.005


Ax = -0.353
Bx = 0.353
Ay = 0.3
By = 0.3
C = 3*pi/8
m = 100
p = 2000
var_list = [-0.16722450068844916, 0.16722450106597264, 0.20752089366702867, 2.032682657322178, 2.032682661716591]
g = 10
alpha1 = 3*pi/2 - var_list[3]
alpha2 = 3*pi/2


def F(var_list):
    f = [0] * 5
    f[0] = var_list[0] + var_list[2] * cos(3 * pi / 2 - var_list[3]) - Ax
    f[1] = var_list[1] + var_list[2] * cos(3 * pi / 2 + var_list[4]) - Bx
    f[2] = var_list[2] + var_list[2] * sin(3 * pi / 2 - var_list[3]) - Ay
    f[3] = (var_list[3] + var_list[4]) * var_list[2] + (var_list[1] - var_list[0]) - C
    f[4] = var_list[2] + var_list[2] * sin(3 * pi / 2 + var_list[4]) - Ay
    return f


def draw(ax, var_list):
    line2, = ax.plot([Ax, Bx], [Ay, By], color='black')
    line4, = ax.plot([-1, 1], [0, 0], color='red')
    line4.set_dashes([2, 2, 4, 2])
    line4.set_dash_capstyle("round")
    line3, = ax.plot([var_list[0], var_list[1]], [0, 0], color='black')
    ellipse_1 = matplotlib.patches.Arc((var_list[0], var_list[2]), 2 * var_list[2], 2 * var_list[2], angle=0,
                                       theta1=degrees(3 * pi / 2 - var_list[3]),
                                       theta2=degrees(3 * pi / 2), color='black')
    ellipse_2 = matplotlib.patches.Arc((var_list[1], var_list[2]), 2 * var_list[2], 2 * var_list[2], angle=-90,
                                       theta1=0,
                                       theta2=degrees(var_list[4]), color='black')
    ax.add_patch(ellipse_1)
    ax.add_patch(ellipse_2)
    camera.snap()


delta_t = 0
Vy = 0

x = np.linspace(-2, 2, 100)

fig, ax = plt.subplots()
camera = Camera(fig)

ax.set(xlim=(-0.5, 0.5), ylim=(-0.5, 0.5))


while delta_t < TIMECONST:
    delta_t += STEP1
    Vy += (1 / m) * (p*(var_list[1] - var_list[0])-m*g) * STEP1
    Ay += Vy * STEP1
    By = Ay
    while np.linalg.norm(F(var_list)) > 0.0001:
        f = F(var_list)
        for j in range(5):
            var_list[j] -= STEP2 * f[j]

    draw(ax, var_list)

animation = camera.animate(interval=10)
plt.show()
