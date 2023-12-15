from math import sqrt


class Vect:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return f"({self.x}, {self.y}, {self.z})"

    def __add__(self, other):
        x = self.x + other.x
        y = self.y + other.y
        z = self.z + other.z
        return Vect(x, y, z)

    def __sub__(self, other):
        x = self.x - other.x
        y = self.y - other.y
        z = self.z - other.z
        return Vect(x, y, z)

    def __mul__(self, other):
        x = self.x * other
        y = self.y * other
        z = self.z * other
        return Vect(x, y, z)

    def __truediv__(self, other):
        x = self.x / other
        y = self.y / other
        z = self.z / other
        return Vect(x, y, z)

    def __neg__(self):
        return self * (-1)

    def __iadd__(self, other):
        new_self = self + other
        self.x = new_self.x
        self.y = new_self.y
        self.z = new_self.z
        return self

    def __isub__(self, other):
        new_self = self - other
        self.x = new_self.x
        self.y = new_self.y
        self.z = new_self.z
        return self

    def __imul__(self, other):
        new_self = other * self
        self.x = new_self.x
        self.y = new_self.y
        self.z = new_self.z
        return self

    def __itruediv__(self, other):
        new_self = self / other
        self.x = new_self.x
        self.y = new_self.y
        self.z = new_self.z
        return self

    def copy(self):
        return Vect(self.x, self.y, self.z)

    def norm_square(self):
        return self.x**2 + self.y**2 + self.z**2

    def norm(self):
        return sqrt(self.norm_square())

    def normalised(self):
        return self/self.norm()

    def do_normalise(self):
        self.__itruediv__(self.norm())
        return self


def dot(v, w):
    return v.x*w.x + v.y*w.y + v.z*w.z


def cross(v, w):
    x = v.y*w.z - v.z*w.y
    y = v.z*w.x - v.x*w.z
    z = v.x*w.y - v.y*w.x
    return Vect(x, y, z)


class RGBColor:
    def __init__(self, r, g, b):
        self.r = r
        self.g = g
        self.b = b

    def __str__(self):
        return f"({self.r}, {self.g}, {self.b})"

    def __add__(self, other):
        r = self.r + other.r
        g = self.g + other.g
        b = self.b + other.b
        return RGBColor(r, g, b)

    def __mul__(self, other):
        r = self.r * other
        g = self.g * other
        b = self.b * other
        return RGBColor(r, g, b)

    def __truediv__(self, other):
        r = self.r / other
        g = self.g / other
        b = self.b / other
        return RGBColor(r, g, b)

    def __iadd__(self, other):
        new_self = self + other
        self.r = new_self.r
        self.g = new_self.g
        self.b = new_self.b
        return self

    def __imul__(self, other):
        new_self = self * other
        self.r = new_self.r
        self.g = new_self.g
        self.b = new_self.b
        return self

    def __itruediv__(self, other):
        new_self = self / other
        self.r = new_self.r
        self.g = new_self.g
        self.b = new_self.b
        return self

    def rounded(self):
        r = int(self.r//1)
        g = int(self.g//1)
        b = int(self.b//1)
        return RGBColor(r, g, b)

    def do_round(self):
        new_self = self.rounded()
        self.r = new_self.r
        self.g = new_self.g
        self.b = new_self.b
        return self

    def write_color(self):
        return str(self.rounded()) + '\n'


class Ray:
    def __init__(self, origin, direction):
        self.origin = origin
        self.direction = direction

    def __str__(self):
        return f"{self.origin}, {self.direction}"

    def copy(self):
        origin = self.origin.copy()
        direction = self.direction.copy()
        return Ray(origin, direction)

    def point_at_t(self, t):
        return self.origin + self.direction * t
