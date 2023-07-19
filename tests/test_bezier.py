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
              [3.0, 1.0], [4.0, 1.0],
              [5.0, 0.5], [6.0, -1.0]]

    bezier = BaseBezier(points = points)
    split_curve_point = []
    split_t = 0.3

    # Split curve
    bezier_split = bezier.de_casteljau(t=split_t)
    cp1 = np.array(bezier_split[0].points)
    cp2 = np.array(bezier_split[1].points)

    tangent = bezier.tangent(t=split_t)
    point = bezier.get_point(t=split_t)
    derivative = bezier.derivative(split_t)
    t = bezier.get_t(point=point)

    coordinates = np.array(bezier.get_coordinates(500))
    l = []
    for i, pi in enumerate(coordinates[1:]):
        l.append(sum([(p1 - p2) ** 2 for p1, p2 in zip(coordinates[i], pi)]) ** 0.5)

    print('Check length')
    print(f'{sum(l)=}, {bezier.get_length()=}')

    # Bezier curve passing through points
    new_points = bezier.get_coordinates(20)
    thr_pnts = BezierThroughPoints(points=new_points)

    ts = linspace(0, 1, 20)

    thr_pnts_crd = []
    for b in thr_pnts.curves:
        thr_pnts_crd.append(b.get_points(ts))

    bpoints = np.array(bezier.get_points(ts))
    points1 = np.array(bezier_split[0].get_points(ts))
    points2 = np.array(bezier_split[1].get_points(ts))
    points = np.array(points)
    new_points = np.array(new_points)
    thr_pnts_crd = np.array(thr_pnts_crd)
    thr_pnts_crd = thr_pnts_crd.reshape((380, 2))

    fig, ax = plt.subplots()
    ax.plot(bpoints[:,0], bpoints[:,1], label='Initial Curve', marker='', color='k')
    ax.plot(points1[:,0], points1[:,1], label='First split curve', marker='', color='r')
    ax.plot(points2[:,0], points2[:,1], label='Second split curve', marker='', color='b')
    ax.scatter(points[:,0], points[:,1], label='Control Points')
    ax.scatter(cp1[:,0], cp1[:,1], label='1st curve CP', marker='o', color='r')
    ax.scatter(cp2[:,0], cp2[:,1], label='2nd curve CP', marker='o', color='b')
    ax.plot(thr_pnts_crd[:,0], thr_pnts_crd[:,1], label='Throught Points', marker='', color='g')
    ax.scatter(new_points[:, 0], new_points[:, 1], color='k', marker='o')

    for i, p in enumerate(points[1:]):
        ax.plot([points[i, 0], p[0]], [points[i, 1], p[1]], marker='', color='k', linestyle='--')

    ax.grid(True)
    fig.legend(loc='upper right')
    plt.show()
