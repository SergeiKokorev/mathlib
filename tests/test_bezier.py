import os
import sys
import matplotlib.pyplot as plt
import numpy as np


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from misc.misc import linspace
from bezier.bezier import BezierThroughPoints, BaseBezier


if __name__ == "__main__":

    points = [[0.0, 0.0], [2.0, 2.0], 
              [3.0, 1.0], [4.0, 1.0]]

    bezier = BaseBezier(points = points)
    split_curve_point = []
    split_t = 0.3

    # Split curve
    bezier_split = bezier.de_casteljau(t=split_t)
    cp1 = np.array(bezier_split[0].points)
    cp2 = np.array(bezier_split[1].points)

    tangent = bezier.tangent(t=split_t)
    point = bezier.get_point(t=split_t)

    ts = linspace(0, 1, 100)
    bpoints = np.array(bezier.get_points(ts))
    points1 = np.array(bezier_split[0].get_points(ts))
    points2 = np.array(bezier_split[1].get_points(ts))
    points = np.array(points)

    fig, ax = plt.subplots()
    ax.plot(bpoints[:,0], bpoints[:,1], label='Initial Curve', marker='', color='k')
    ax.plot(points1[:,0], points1[:,1], label='First split curve', marker='', color='r')
    ax.plot(points2[:,0], points2[:,1], label='Second split curve', marker='', color='b')
    ax.scatter(points[:,0], points[:,1], label='Control Points')
    ax.scatter(cp1[:,0], cp1[:,1], label='1st curve CP', marker='o', color='r')
    ax.scatter(cp2[:,0], cp2[:,1], label='2nd curve CP', marker='o', color='b')
    ax.plot([points[0,0], points[1,0]], [points[0,1], points[1,1]], marker='', color='k')
    ax.plot([points[1,0], points[2,0]], [points[1,1], points[2,1]], marker='', color='k')
    ax.plot([points[-1,0], points[-2,0]], [points[-1,1], points[-2,1]], marker='', color='k')

    ax.grid(True)
    fig.legend(loc='upper right')
    plt.show()
