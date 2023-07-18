from __future__ import annotations

import os
import sys


from copy import deepcopy
from typing import List, Tuple


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))

sys.path.append(DIR)

from misc.misc import *
from matrix.matrix import *
from consts.gaussian_quadrature import GQ as gq


class BaseBezier:
    """
    Class performs basic operation on Bezier object
    Properties
    ----------
    points: list of Bezier curve points: list
    degree: Bezier curve degree: int = number 0f control points - 1
    dimention: Bezier curve dimention: int = number of point's coordinates
    num_point: number of interpolation points of Bezier curve: int (default 50)
    first_point: First point of Bezier curve: list
    last_point: Last point of Bezier curve: list

    Methods
    -------
    get_polynomial_coefficients() -> list
        Returns the polynomial coefficients of the Bezier curve
        list[a0, a1, a2, ..., an], where a0 is a free term, a coefficient at x**n
    get_point(t: float=0.0) -> list
        Returns point at t. list [x, y, z]
    get_points(t: list) -> list
        Returns a list of points at ts. list[[x0, y0, z0], [x1, y1, z1], ...]
    get_length(n: int=5, t1: float=0.0, t2: float=1.0) -> float
        Returns a length of Bezier curve between points t1 and t2, 
        computed by using Gaussian Quadrature, n is number of points 
        used for computing length
    get_coordinates(npoints: int=0, start: float=0.0, stop: float=1.0) -> list
        Returns Cartesian coordinates of points lying between start and stop
        npoints is number of points
    get_t(point: tuple([0, 0.0]), start: float=0.0, stop=1.0, eps: float=10**-3) -> float
        Returns Bezier curve parameter t at point point. Where point is a tuple.
        point[0] is a number of coordinate (ex. point[0]=0 is first coordinate X)
        point[1] this is a coordinate itselfin which t will be computed
        start and stop are the interval over which the t will be computed
        eps is a computation accurance
    def tangetn(t: float=0.0) -> list
        Returns normalized tangent vector to the curve at the point t
    get_length_point(length: float=0.0, t1: float=0.0, t2: float=1.0)
        Returns point at particular Bezier length
        length is the length at this point will be computed
        t1 and t2 are the interval of computation
    def de_casteljau(t: float = 0.5) -> List[BaseBeizer]
        Split Bezier at point t and returns two new BaseBezier curve
        t is a parameter of Bezier,
        return List[BaseBezier]
    def __de_casteljau(t: float = 0.5) -> Tuple[List]
        split curve and return tuple of new control points of two curves
        return tuple(b1, b2)
    """

    def __init__(self, points: list, *args, **kwargs):
        self.__points = points
        self.__degree = len(self.__points) - 1
        self.__space_dimention = len(self.__points[0])
        self.__npoints = kwargs.get('npoints', 50)
        self.__first_point = self.__points[0]
        self.__last_point = self.__points[-1]

    @property
    def points(self):
        return self.__points

    @property
    def degree(self):
        return self.__degree

    @property
    def dimention(self):
        return self.__space_dimention

    @property
    def num_points(self):
        return self.__npoints

    @property
    def first_point(self):
        return self.__first_point

    @property
    def last_point(self):
        return self.__last_point

    def _binomial_coefficient(self, n: int, i: int) -> float:
        return factorial(n) / (factorial(i) * factorial(n - i))

    def _Bernstein_polynomial(self, n: int, i: int, t: float) -> float:
        return self._binomial_coefficient(n, i) * t ** i * (1 - t) ** (n - i)

    def get_polynomial_coefficients(self) -> list:

        c = []
        n = self.degree
        cols = self.dimention
        for j in range(0, n + 1):
            ci = [0] * cols
            binom = factorial(n) / factorial(n - j)
            for i in range(j + 1):
                polynom = ((-1) ** (i + j)) / (factorial(i) * factorial(j - i))
                for col in range(cols):
                    ci[col] += binom * polynom * self.points[i][col]
            c.append(ci)
        return c
                    
    def get_point(self, t: float=0) -> list:
        return self.__de_casteljau(t)[0][-1]

    def get_points(self, t: list) -> list:
        return [self.get_point(ti) for ti in t]

    def get_length(self, t1: float=0.0, t2: float=1.0, n: int=5):
        w = gq[n][0]
        x = gq[n][1]
        length = 0.0
        for i in range(n):
            ti = (t2 - t1) * x[i] / 2 + (t2 - t1) / 2
            ft = sum([p ** 2 for p in self.derivative(t=ti)]) ** 0.5
            length += w[i] * ft
        length *= (t2 - t1) / 2
        return length

    def get_coordinates(self, npoints: int = 10, start: float = 0.0, stop: float = 1.0) -> list:
        points = []
        t = linspace(start, stop, npoints)
        for ti in t:
            points.append(self.__de_casteljau(ti)[0])
        
        return points

    def get_t(self, point: list, t0: float = 0.0, t1: float = 1.0, eps: float = 10 ** -4, it: int = 0):

        max_it = 100
        it = it

        left, right = self.get_point(t0), self.get_point(t1)
        t = linspace(t0, t1, 500)
        mid = self.get_point(t[int(len(t) / 2)])

        error = sum([(mi - pi) ** 2 for mi, pi in zip(mid, point)]) ** 0.5
        lerror = sum([(mi - pi) ** 2 for mi, pi in zip(left, point)]) ** 0.5
        rerror = sum([(mi - pi) ** 2 for mi, pi in zip(right, point)]) ** 0.5

        if error <= eps : return t[int(len(t) / 2)]
        if lerror <= eps : return t0
        if rerror <= eps : return t1
        if it >= max_it : return False

        if lerror > rerror:
            t0 = t[int(len(t) / 2)]
            t1 = t[len(t) - 1]
        else:
            t0 = t[0]
            t1 = t[int(len(t) / 2)]
        return self.get_t(point, t0, t1, eps, it+1)

    def get_length_point(self, length: float=0.0, t1: float=0.0, t2: float=1.0):
        eps = 10 ** -3

        if (length - 0) <= eps:
            return 0.0, self.get_point(t=0.0)
        elif abs(length - self.get_length()) <= eps:
            return 1.0, self.get_point(t=1.0)

        cur_len = 0.0
        count = 0.0
        max_count = 10 ** 5

        while abs(cur_len - length) > eps:
            if length > cur_len:
                t2 = t2 + 0.1 * t2
            elif length < cur_len:
                t2 = t2 - 0.1 * t2

            cur_len = self.get_length(t1=t1, t2=t2)
            count += 1
            if count > max_count:
                return False
        return t2, self.get_point(t=t2)

    def de_casteljau(self, t: float = 0.5) -> List[BaseBezier]:

        beta = [[c for c in self.points]]
        b1 = []
        b2 = []
        n = self.degree + 1

        for i in range(1, n):
            beta.insert(i, [[0.0] * self.dimention] * n)

        for j in range(1, n):
            for i in range(0, n - j):
                beta[j][i] = [b1 * (1 - t) + b2 * t for b1, b2 in zip(beta[j-1][i], beta[j-1][i+1])]

        for i, b in enumerate(beta):
            for j, be in enumerate(b):
                if j == 0:
                    b1.append(be)
                if (i + j) == n - 1:
                    b2.append(be)

        return [BaseBezier(points=b1[::1]), BaseBezier(points=b2[::-1])]

    def tangent(self, t: float = 0) -> list:
        points = self.__de_casteljau(t)
        tan = [p1 - p2 for p1, p2 in zip(points[0][-2], points[1][1])]
        magnitude = sum([p ** 2 for p in tan]) ** 0.5
        return [p / magnitude for p in tan]

    def __de_casteljau(self, t: float = 0.5) -> Tuple[List]:

        beta = [[c for c in self.points]]
        b1 = []
        b2 = []
        n = self.degree + 1
        for i in range(1, n):
            beta.insert(i, [[0.0] * self.dimention] * n)
        for j in range(1, n):
            for i in range(0, n - j):
                beta[j][i] = [b1 * (1 - t) + b2 * t for b1, b2 in zip(beta[j-1][i], beta[j-1][i+1])]
        for i, b in enumerate(beta):
            for j, be in enumerate(b):
                if j == 0:
                    b1.append(be)
                if (i + j) == n - 1:
                    b2.append(be)
        return b1, b2[::-1]

    def __repr__(self):
        return 'BaseBezier\ndegree {degree};\n' \
               'Points: {points}\n'.format(degree=self.degree, points=self.points)


class BezierThroughPoints(BaseBezier):
    """
    Class creates piecewise Bzeier curve through points and provides basic operations
    Properties
    ----------
    first_point: First point of the Bezier curve: list
    last_point: Last point of the Bezier curve: list
    degree: Bezier curve degree: int=3
    lenght: Summary length of all pieces of the Bezier curve: float
    lengths: List of all pieces length: list
    A and B: A equation coefficient matrices to compute Bezier curve control points: list
    num_curves: Number of pieces of curve

    Methods
    -------
    interpolate(points: tuple([0, list])) -> list
        Return interpolated coordinates of the Bezier curve
        points is a tuple. points[0] a coordinate index
        points[1] is a list of coordinates
    norm_length_point(norm_length: float=0.0) -> list
        Returns a coordinate of point at normed length = norm_length
        [x, y, z]
    get_coordinates(npoints: int=0) -> list
        Returns coordinates of points. npoints is a number of points,
        if is None then the self.num_points(default=50) is taken
        [[x0, y0, z0], [x1, y1, z1], ...]
    derivative(norm_length: float=0.0) -> list
        Returns a derivative of the Bezier curve at the point 
        corresponding the normalized length=norm_length
    derivatives(norm_length: list=[0.0, 1.0]) -> list
        Returns derivatives at the points corresponding the normalized
        lengths in the list norm_length
    """

    def __init__(self, points: list, *args, **kwargs):
        super().__init__(points=points, *args, **kwargs)
        self.curves: List[BaseBezier] = []
        self.__degree = kwargs.get('degree', 3)
        self.__dim = 2 + (len(self.points) - 2) * 2
        self.__A = [[0] * self.__dim for d in range(self.__dim)]
        self.__B = [[0] * self.dimention for d in range(self.__dim)]
        self.__number_of_curves = len(self.points) - 1
        for points in self.__control_points():
            self.curves.append(BaseBezier(points=points))

        self.__first_point = self.curves[0].points[0]
        self.__last_point = self.curves[-1].points[-1]
        self.__curve_lengths = [curve.get_length() for curve in self.curves]
        self.__length = sum(self.__curve_lengths)
        self._control_points = self.__control_points()

    @property
    def first_point(self):
        return self.__first_point

    @property
    def last_point(self):
        return self.__last_point

    @property
    def degree(self):
        return self.__degree

    @property
    def length(self):
        return self.__length

    @property
    def A(self) -> list:

        self.__A[0][:2] = [2, -1]
        self.__A[self.__dim - 1][self.__dim - 2:] = [-1, 2]

        j = 0
        for i in range(self.__dim):
            if j < self.__dim - 2:
                self.__A[i + 1][j + 1 : j + 3] = [1, 1]
                self.__A[i + int(self.__dim / 2)][j : j + 4] = [1, -2, 2, -1]
            j += 2

        return self.__A

    @property
    def B(self) -> list:

        self.__B[0] = self.points[0]
        self.__B[self.__dim - 1] = self.points[len(self.points) - 1]
        self.__B[1:len(self.points) - 1] = [
            [2 * p[i] for i in range(self.dimention)] for p in self.points[1:len(self.points) - 1]
            ]

        return self.__B

    @property
    def num_curves(self):
        return self.__number_of_curves

    @property
    def lengths(self):
        return self.__curve_lengths

    @property
    def control_points(self):
        return self._control_points

    def get_point(self, t: float) -> list:

        pass

    def get_coordinates(self, npoints: int = 10) -> list:

        pass

    def derivative(self, t: float = 0.5) -> list:
        
        pass

    def derivatives(self, t: list = [0.0, 1.0]) -> list:
        
        pass

    def tangent_line(self, t: float = 0.5):

        pass

    def __control_points(self) -> list:

        print('Compute control points')

        cp = matmul(inverse(self.A), self.B)
        control_points = []

        for i in range(self.num_curves):
            tmp = []
            tmp.append(self.points[i])
            for p in cp[2 * i:2 * i + 2]:
                tmp.append(p)
            tmp.append(self.points[i + 1])
            control_points.append(tmp)

        return control_points
