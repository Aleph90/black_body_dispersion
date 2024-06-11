from math import sqrt


# Geometry and linear algebra
class Vect3D:
    def __init__(self, x: float, y: float, z: float):
        self._entries = [x, y, z]

    @classmethod
    def cast(cls, v: 'Vect3D') -> 'Vect3D':
        return cls(*v.entries)

    @classmethod
    def null(cls):
        return cls(0, 0, 0)

    def __eq__(self, other) -> bool:
        return self._entries == other.entries

    def __copy__(self) -> 'Vect3D':
        return Vect3D(*self._entries)

    def copy(self) -> 'Vect3D':
        return self.__copy__()

    # Representations and type conversion
    def __repr__(self) -> str:
        return f'Vect3D{tuple(self._entries)}'

    def __str__(self) -> str:
        return f'{tuple(self._entries)}'

    def __bool__(self) -> bool:
        return bool([e for e in self._entries if e])

    # Access
    @property
    def entries(self) -> list[float]:
        return self._entries

    @property
    def x(self) -> float:
        return self._entries[0]

    @property
    def y(self) -> float:
        return self._entries[1]

    @property
    def z(self) -> float:
        return self._entries[2]

    def __getitem__(self, key: int) -> float:
        return self._entries[key]

    @x.setter
    def x(self, value: float) -> None:
        self._entries[0] = value

    @y.setter
    def y(self, value: float) -> None:
        self._entries[1] = value

    @z.setter
    def z(self, value: float) -> None:
        self._entries[2] = value

    def __setitem__(self, key: int, value: float) -> None:
        self._entries[key] = value

    # Linear algebra
    def __add__(self, other: 'Vect3D') -> 'Vect3D':
        return Vect3D(*[self._entries[i] + other[i] for i in range(3)])

    def __sub__(self, other: 'Vect3D') -> 'Vect3D':
        return self + (-other)

    def __rmul__(self, scalar: float) -> 'Vect3D':
        return Vect3D(*[scalar * entry for entry in self._entries])

    def __mul__(self, scalar: float) -> 'Vect3D':
        return scalar * self

    def __truediv__(self, scalar: float) -> 'Vect3D':
        return (1/scalar) * self

    def __neg__(self) -> 'Vect3D':
        return (-1) * self

    # Euclidean structure
    def dot_with(self, other: 'Vect3D') -> float:
        return sum([self._entries[i] * other[i] for i in range(3)])

    def norm_squared(self) -> float:
        return self.dot_with(self)

    def norm(self) -> float:
        return sqrt(self.norm_squared())

    def normalised(self) -> 'Vect3D':
        return (1/self.norm()) * self

    def cross_with(self, other: 'Vect3D') -> 'Vect3D':
        return Vect3D(self.y * other.z - self.z * other.y,
                      self.z * other.x - self.x * other.z,
                      self.x * other.y - self.y * other.x)


class Matrix3by3:
    def __init__(self, c1: Vect3D | list[int], c2: Vect3D | list[int], c3: Vect3D | list[int]):
        if isinstance(c1, list):
            c1 = Vect3D(*c1)
        if isinstance(c2, list):
            c2 = Vect3D(*c2)
        if isinstance(c3, list):
            c3 = Vect3D(*c3)
        self._cols = [c1.copy(), c2.copy(), c3.copy()]

    def __eq__(self, other: 'Matrix3by3') -> bool:
        return self._cols == other._cols

    def __copy__(self) -> 'Matrix3by3':
        return Matrix3by3(*[c for c in self._cols])

    def __deepcopy__(self) -> 'Matrix3by3':
        return Matrix3by3(*[c.copy() for c in self._cols])

    # Representation and type conversion
    def __repr__(self) -> str:
        return f'Matrix3by3({", ".join([c.__repr__() for c in self._cols])})'

    def __str__(self) -> str:
        return f'[{[c.x for c in self._cols]},\n{[c.y for c in self._cols]},\n{[c.z for c in self._cols]}]'

    def __bool__(self) -> bool:
        return bool([c for c in self._cols if c])

    # Access
    def __getitem__(self, pos: tuple[int, int]) -> float:
        i, j = pos
        return self._cols[j][i]

    def __setitem__(self, pos: tuple[int, int], val: float) -> None:
        i, j = pos
        self._cols[j][i] = val

    # Algebra
    def __rmul__(self, scalar: float) -> 'Matrix3by3':
        return Matrix3by3(*[scalar * col for col in self._cols])

    def __mul__(self, other):
        if isinstance(other, float):
            product = other * self
        elif isinstance(other, Vect3D):
            product = sum([other[i]*self._cols[i] for i in range(3)], start=self._cols[0].null())
        elif isinstance(other, Matrix3by3):
            product = Matrix3by3(*[self * col for col in other._cols])
        else:
            raise TypeError('What are you trying to multiply this matrix by?')
        return product

    def det(self) -> float:
        return dot(self._cols[0], cross(self._cols[1], self._cols[2]))


def dot(v: Vect3D, w: Vect3D) -> float:
    return v.dot_with(w)


def norm_squared(v: Vect3D) -> float:
    return v.norm_squared()


def norm(v: Vect3D) -> float:
    return v.norm()


def normalised(v: Vect3D) -> Vect3D:
    return v.normalised()


def cross(v: Vect3D, w: Vect3D) -> Vect3D:
    return v.cross_with(w)


def det(m: Matrix3by3) -> float:
    return m.det()


def is_basis(u: Vect3D, v: Vect3D, w: Vect3D, precision: int = 3):
    return abs(det(Matrix3by3(u, v, w))) > 10**(-precision)


def is_positive_frame(u: Vect3D, v: Vect3D, w: Vect3D):
    return det(Matrix3by3(u, v, w)) > 0


class Ray:
    def __init__(self, origin: Vect3D, direction: Vect3D):
        self.origin = origin.copy()
        self.direction = direction.copy()

    def __eq__(self, other: 'Ray') -> bool:
        return self.origin == other.origin and self.direction == other.direction

    def __copy__(self) -> 'Ray':
        return Ray(self.origin, self.direction)

    def copy(self) -> 'Ray':
        return self.__copy__()

    def __repr__(self) -> str:
        return f'Ray({str(self.origin)}, {str(self.direction)})'

    def __str__(self) -> str:
        return f'({str(self.origin)}, {str(self.direction)})'

    def point_at(self, t: float) -> Vect3D:
        return self.origin + t * self.direction


# Color
class RGBColor(Vect3D):
    def __init__(self, r: float, g: float, b: float):
        super().__init__(r, g, b)

    def __eq__(self, other) -> bool:
        return isinstance(other, RGBColor) and super().__eq__(other)

    def __copy__(self) -> 'RGBColor':
        return RGBColor(*self._entries)

    def copy(self) -> 'RGBColor':
        return self.__copy__()

    # Representations and type conversion
    def __repr__(self) -> str:
        return f'RGBColor{tuple(self._entries)}'

    def __str__(self) -> str:
        return super().__str__()

    def __bool__(self) -> bool:
        return super().__bool__()

    # Access
    @property
    def entries(self) -> list[float]:
        return super().entries

    @property
    def r(self) -> float:
        return self.x

    @property
    def g(self) -> float:
        return self.y

    @property
    def b(self) -> float:
        return self.z

    def __getitem__(self, key: int) -> float:
        return super().__getitem__(key)

    @r.setter
    def r(self, value: float) -> None:
        self.x = value

    @g.setter
    def g(self, value: float) -> None:
        self.y = value

    @b.setter
    def b(self, value: float) -> None:
        self.z = value

    def __setitem__(self, key: int, value: float) -> None:
        super().__setitem__(key, value)

    # Algebra
    def __add__(self, other: 'RGBColor') -> 'RGBColor':
        return self.cast(super().__add__(other))

    def __sub__(self, other: 'RGBColor') -> 'RGBColor':
        return self.cast(super().__sub__(other))

    def __rmul__(self, scalar: float) -> 'RGBColor':
        return self.cast(super().__rmul__(scalar))

    def __mul__(self, scalar: float) -> 'RGBColor':
        return self.cast(super().__mul__(scalar))

    def __truediv__(self, scalar: float) -> 'RGBColor':
        return self.cast(super().__truediv__(scalar))

    def __neg__(self) -> 'RGBColor':
        return self.cast(super().__neg__())


class SumFun:
    def __init__(self, *components):
        self._components = list(components)

    # Representation and type conversion
    def __repr__(self):
        return f'SumFun{tuple(self._components)}'

    def __str__(self):
        return ' + '.join([str(comp) for comp in self._components])

    # Function call
    def __call__(self, *args) -> float:
        return sum([component(*args) for component in self._components])

    # Algebra
    def __add__(self, other: 'SumFun') -> 'SumFun':
        return SumFun(*self._components + other._components)


class MathFunctionWithLookup:
    def __init__(self, fun):
        self._fun = fun
        self._lookup = {}

    def __call__(self, *args):
        if args not in self._lookup:
            self._lookup[args] = self._fun(*args)
        return self._lookup[args]
