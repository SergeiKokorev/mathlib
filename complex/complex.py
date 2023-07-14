from __future__ import annotations


from typing import Any
from dataclasses import dataclass
from copy import deepcopy


@dataclass
class Vector:
    x: float | int = 0.0
    y: float | int = 0.0
    z: float | int = 0.0

    @property
    def norm(self) -> float:
        return self.x * self.x + self.y * self.y + self.z * self.z
    
    @property
    def magnitude(self) -> float:
        return self.norm ** 0.5
    
    def cross(self, var) -> float:
        if not isinstance(var, (Vector, float, int)):
            raise TypeError(f'unsupported operand type(s) for cross product: {type(self)} and {type(var)}')
        if isinstance(var, (float, int)):
            var = Vector(var, 0.0, 0.0)
        return Vector(
            self.y * var.z - self.z * var.y,
            self.z * var.x - self.x * var.z,
            self.x * var.y - self.y * var.x
        )
    
    def normalize(self) -> None:
        m = self.magnitude
        self.x /= m
        self.y /= m
        self.z /= m

    def __add__(self, var) -> Vector:
        if not isinstance(var, (int, float, Vector)):
            raise TypeError(f'unsupported operand type(s) for +: {type(self)} and {type(var)}')
        if isinstance(var, (int, float)):
            var = Vector(var, var, var)
        return Vector(self.x + var.x, self.y + var.y, self.z + var.z)
    
    def __radd__(self, var) -> Vector:
        return self.__add__(var)
    
    def __iadd__(self, var) -> Vector:
        return self.__add__(var)

    def __sub__(self, var) -> Vector:
        if not isinstance(var, (int, float, Vector)):
            raise TypeError(f'unsupported operand type(s) for -: {type(self)} and {type(var)}')
        if isinstance(var, (int, float)):
            var = Vector(var, var, var)
        return Vector(self.x - var.x, self.y - var.y, self.z - var.z)
    
    def __rsub__(self, var) -> Vector:
        if not isinstance(var, (int, float, Vector)):
            raise TypeError(f'unsupported operand type(s) for -: {type(self)} and {type(var)}')
        if isinstance(var, (int, float)):
            var = Vector(var, var, var)
        return Vector(var.x - self.x, var.y - self.y, var.z - self.z)
    
    def __isub__(self, var) -> Complex:
        return self.__sub__(var)

    def __mul__(self, var) -> float | Vector:
        if not isinstance(var, (int, float, Vector)):
            raise TypeError(f'unsupported operand type(s) for dot product: {type(self)} and {type(var)}')
        if isinstance(var, (int, float)):
            return Vector(self.x * var, self.y * var, self.z * var)
        return self.x * var.x + self.y * var.y + self.z * var.z

    def __rmul__(self, var) -> Vector:
        return self.__mul__(var)
    
    def __imul__(self, var) -> Vector:
        return self.__mul__(var)

    def __truediv__(self, var) -> Vector:

        if not isinstance(var, (float, int)):
            raise TypeError(f'unsupported operand type(s) for /: {type(self)} and {type(var)}')
        return Vector(self.x / var, self.y / var, self.z / var)

    def __str__(self) -> str:
        return f'{self.x}i + {self.y}j + {self.z}k'


@dataclass
class Complex: 
    Re: int | float = 1.0
    Im: int | float = 0.0

    @property
    def magnitude(self):
        return (self.Re * self.Re + self.Im * self.Im) ** 0.5

    def conjugate(self) -> None:
        self.Im = (-1) * self.Im

    def inverse(self) -> None:
        return Complex(self.Re / self.magnitude, self.Im / self.magnitude)

    def __add__(self, var) -> Complex:
        if not isinstance(var, (float, int, Complex)): 
            raise TypeError(f'unsupported operand type(s) for +: {type(self)} and {type(var)}')
        elif isinstance(var, (float, int)):
            var = Complex(var, 0.0)
        return Complex(self.Re + var.Re, self.Im + var.Im)
    
    def __radd__(self, var) -> Complex:
        return self.__add__(var)

    def __iadd__(self, var) -> Complex:
        return self.__add__(var)

    def __sub__(self, var) -> Complex:
        if not isinstance(var, (float, int, Complex)): 
            raise TypeError(f'unsupported operand type(s) for -: {type(self)} and {type(var)}')
        elif isinstance(var, (float, int)):
            var = Complex(var, 0.0)
        return Complex(self.Re - var.Re, self.Im - var.Im)
    
    def __rsub__(self, var) -> Complex:
        if not isinstance(var, (float, int, Complex)): 
            raise TypeError(f'unsupported operand type(s) for -: {type(self)} and {type(var)}')
        elif isinstance(var, (float, int)):
            var = Complex(var, 0.0)
        return Complex(var.Re - self.Re, var.Im - self.Im)

    def __isub__(self, var) -> Complex:
        return self.__sub__(var)

    def __mul__(self, var) -> Complex:
        if not isinstance(var, (float, int, Complex)): 
            raise TypeError(f'unsupported operand type(s) for *: {type(self)} and {type(var)}')
        elif isinstance(var, (int, float)): 
            var = Complex(var, 0.0)
        return Complex(self.Re * var.Re - self.Im * var.Im, self.Re * var.Im + self.Im * var.Re)
    
    def __rmul__(self, var) -> Complex:
        if not isinstance(var, (float, int, Complex)): 
            raise TypeError(f'unsupported operand type(s) for *: {type(self)} and {type(var)}')
        elif isinstance(var, (int, float)): 
            var = Complex(var, 0.0)
        return Complex(var.Re * self.Re - var.Im * self.Im, var.Re * self.Im + var.Im * self.Re)
    
    def __imul__(self, var) -> Complex:
        return self.__mul__(var)
    
    def __itruedive__(self, var) -> Complex:
        return self.__truediv__(var)

    def __truediv__(self, var) -> Complex:
        if not isinstance(var, (float, int, Complex)): 
            raise TypeError(f'unsupported operand type(s) for /: {type(self)} and {type(var)}')
        elif isinstance(var, (int, float)): 
            var = Complex(var, 0.0)
        return Complex(
            (var.Re * self.Re + var.Im * self.Im) / self.magnitude,
            (var.Im * self.Re - var.Re * self.Im) / self.magnitude
        )
    
    def __rtruediv__(self, var) -> Complex:
        if not isinstance(var, (float, int, Complex)): 
            raise TypeError(f'unsupported operand type(s) for /: {type(self)} and {type(var)}')
        elif isinstance(var, (int, float)): 
            var = Complex(var, 0.0)
        return Complex(
            (self.Re * var.Re + self.Im * var.Im) / var.magnitude,
            (self.Im * var.Re - self.Re * var.Im) / var.magnitude
        )

    def __str__(self) -> str:
        return f'{self.Re} + {self.Im}i'

@dataclass
class DualNumber:
    real: Any = 1.0
    dual: Any = 0.0

    def dnConjugate(self) -> None:
        self.dual = (-1) * self.dual

    @property
    def magnitude(self) -> DualNumber:
        return DualNumber(self.real * self.real, self.real * self.dual - self.dual * self.real)

    def __add__(self, var) -> DualNumber:

        if isinstance(var, (int, float)):
            return DualNumber(self.real + var, self.dual)
        elif isinstance(var, DualNumber):
            return DualNumber(self.real + var.real, self.dual + var.dual)
        else:
            raise TypeError(f'unsupported operand type(s) for +: {type(self)} and {type(var)}')

    def __radd__(self, var) -> DualNumber:
        return self.__add__(var)
    
    def __iadd__(self, var) -> DualNumber:
        return self.__add__(var)

    def __sub__(self, var) -> DualNumber:

        if not isinstance(var, (int, float, DualNumber)):
            raise TypeError(f'unsupported operand type(s) for -: {type(self)} and {type(var)}')

        if isinstance(var, (int, float)):
            var = DualNumber(var, 0.0)

        return DualNumber(self.real - var.real, self.dual - var.dual)

    def __rsub__(self, var) -> DualNumber:
        
        if not isinstance(var, (int, float, DualNumber)):
            raise TypeError(f'unsupported operand type(s) for -: {type(self)} and {type(var)}')

        if isinstance(var, (int, float)):
            var = DualNumber(var, 0.0)

        return DualNumber(var.real - self.real, var.dual - self.dual)

    def __mul__(self, var) -> DualNumber:

        if not isinstance(var, (DualNumber, float, int)):
            raise TypeError(f'unsupported operand type(s) for *: {type(self)} and {type(var)}')

        if isinstance(var, (int, float)):
            var = DualNumber(var, 0.0)

        return DualNumber(self.real * var.real, self.real * var.dual + self.dual * var.real)
            
    def __rmul__(self, var) -> DualNumber:
        if not isinstance(var, (DualNumber, float, int)):
            raise TypeError(f'unsupported operand type(s) for *: {type(self)} and {type(var)}')

        if isinstance(var, (int, float)):
            var = DualNumber(var, 0.0)

        return DualNumber(var.real * self.real, var.real * self.dual + var.dual * self.real)

    def __imul__(self, var) -> DualNumber:
        return self.__mul__(var)

    def __truediv__(self, var) -> DualNumber:

        if not isinstance(var, (int, float, DualNumber)):
            raise TypeError(f'unsupported operand type(s) for /: {type(self)} and {type(var)}')
        elif isinstance(var, (float, int)):
            var = DualNumber(var, 0.0)
        return DualNumber(self.real / var.real, 
                        (self.dual * var.real - self.real * var.dual) / (var.real * var.real))
    
    def __rtruediv__(self, var) -> DualNumber:
        if not isinstance(var, (int, float, DualNumber)):
            raise TypeError(f'unsupported operand type(s) for /: {type(self)} and {type(var)}')
        elif isinstance(var, (float, int)):
            var = DualNumber(var, 0.0)
        return DualNumber(var.real / self.real, 
                        (var.dual * self.real - var.real * self.dual) / (self.real * self.real))

    def __itruediv__(self, var) -> DualNumber:
        return self.__truediv__(var)

    def __str__(self) -> str:
        return f'{self.real} + ({self.dual})e'


@dataclass
class Quaternion(Complex):

    Re: float | int = 1.0
    Im: Vector = Vector(0.0, 0.0, 0.0)

    @property
    def w(self):
        return self.Re
    
    @property
    def x(self):
        return self.Im.x
    
    @property
    def y(self):
        return self.Im.y
    
    @property
    def z(self):
        return self.Im.z

    @w.setter
    def w(self, w: int | float):
        if not isinstance(w, (float, int)):
            raise TypeError(f'Scalar part of an Quaternion must be int or float type (not {type(w)})')
        self.Re = w

    @x.setter
    def x(self, x: int | float):
        if not isinstance(x, (float, int)):
            raise TypeError(f'x component of a vector part of an Quaternion must be int or float type (not {type(x)})')
        self.Im.x = x

    @y.setter
    def y(self, y: int | float):
        if not isinstance(y, (float, int)):
            raise TypeError(f'y component of a vector part of an Quaternion must be int or float type (not {type(y)})')
        self.Im.y = y

    @z.setter
    def z(self, z: int | float):
        if not isinstance(z, (float, int)):
            raise TypeError(f'z component of a vector part of an Quaternion must be int or float type (not {type(z)})')
        self.Im.z = z

    def normalize(self) -> None:
        n = self.magnitude
        self.Re = self.Re / n
        self.Im = self.Im / n

    def inverse(self) -> None:
        m = self.magnitude ** 2
        self.conjugate()
        self.Re = self.Re / m
        self.Im = self.Im / m

    def __add__(self, var) -> Quaternion:
        if not isinstance(var, (Quaternion, int, float)):
            raise TypeError(f'unsupported operand type(s) for +: {type(self)} and {type(var)}')
        
        if isinstance(var, (float, int)):
            return Quaternion(self.Re + var, self.Im + var)
        return Quaternion(var.Re + self.Re, var.Im + self.Im)

    def __radd__(self, var) -> Quaternion:
        return self.__add__(var)

    def __iadd__(self, var) -> Quaternion:
        return self.__add__(var)

    def __sub__(self, var) -> Quaternion:

        if not isinstance(var, (Quaternion, float, int)):
            raise TypeError(f'unsupported operand type(s) for -: {type(self)} and {type(var)}')
        
        if isinstance(var, (float, int)):
            return Quaternion(self.Re - var, self.Im - var)
        
        return Quaternion(self.Re - var.Re, self.Im - var.Im)
    
    def __rsub__(self, var) -> Quaternion:

        if not isinstance(var, (Quaternion, float, int)):
            raise TypeError(f'unsupported operand type(s) for -: {type(self)} and {type(var)}')
        
        if isinstance(var, (float, int)):
            return Quaternion(var - self.Re, var - self.Im)
        
        return Quaternion(var.Re - self.Re, var.Im - self.Im)
    
    def __isub__(self, var) -> Quaternion:
        return self.__sub__(var)

    def __mul__(self, var) -> Quaternion:

        if not isinstance(var, (Quaternion, float, int)):
            raise TypeError(f'unsupported operand type(s) for *: {type(self)} and {type(var)}')
        
        if isinstance(var, (float, int)):
            return Quaternion(self.Re * var, self.Im * var)
        
        return Quaternion(
            self.Re * var.Re - self.Im * var.Im,
            self.Re * var.Im + var.Re * self.Im + self.Im.cross(var.Im)
        )

    def __rmul__(self, var) -> Quaternion:
        return self.__mul__(var)

    def __imul__(self, var) -> Quaternion:
        return self.__mul__(var)

    def __truediv__(self, var) -> Quaternion:
        
        if not isinstance(var, (Quaternion, float, int)):
            raise TypeError(f'unsupported operand type(s) for /: {type(self)} and {type(var)}')
        
        if isinstance(var, (float, int)):
            return Quaternion(self.Re / var, self.Im / var)
        
        q = deepcopy(var)
        q.inverse()
        return Quaternion(
            self.Re * q.Re - self.Im * q.Im,
            self.Re * q.Im + q.Re * self.Im - self.Im.cross(q.Im)
        )

    def __rtruediv__(self, var) -> Quaternion:

        if not isinstance(var, (Quaternion, float, int)):
            raise TypeError(f'unsupported operand type(s) for /: {type(self)} and {type(var)}')
        
        if isinstance(var, (float, int)):
            q = deepcopy(self)
            q.inverse()
            return Quaternion(var * q.Re, var * q.Im)
        
        q = deepcopy(self)
        q.inverse()
        return Quaternion(
            var.Re * q.Re - var.Im * q.Im,
            var.Re * q.Im + q.Re * var.Im - var.Im.cross(q.Im)
        )
    
    def __itruediv__(self, var) -> Quaternion:
        return self.__truediv__(var)

    def __str__(self) -> str:
        return f'{self.Re} + ({self.Im})i'

@dataclass
class DualQuaternion(DualNumber):
    real: Quaternion = Quaternion(1.0, Vector(0.0, 0.0, 0.0))
    dual: Quaternion = Quaternion(0.0, Vector(0.0, 0.0, 0.0))

    @property
    def norm(self):
        q0 = deepcopy(self.real)
        q1 = deepcopy(self.dual)
        q0.conjugate()
        q1.conjugate()

        return DualNumber(
            self.real.magnitude,
            ((q0 * self.dual + q1 * self.real) / (2 * self.real.magnitude)).w
        )
    
    def inverse(self):
        norm = self.norm * self.norm
        q0 = deepcopy(self.real)
        q1 = deepcopy(self.dual)
        q0.conjugate()
        q1.conjugate()
        self.dual = (q1 * norm.real - q0 * norm.dual)
        self.real = q0 / norm.real

    def normalize(self):
        norm = self.norm
        self.dual = (self.dual * norm.real - self.real * norm.dual) / (norm.real * norm.real)
        self.real = self.real / norm.real

    def qConjugate(self):
        self.real.conjugate()
        self.dual.conjugate()

    def conjugate(self):
        self.real.conjugate()
        self.dual.conjugate()
        self.dnConjugate()
