"""This module implements utilities for "hittable" objects (naming convention
taken from the 'Ray tracing in one weekend' tutorial). Every such object
comes with a unit normal vector which represents the *outgoing* direction.

Each class for hittable objects implements methods to detect intersections
with rays, determine the time of intersection (according to the parameterization
of the ray), and record information on the intersection. It also includes
methods to determine whether a point is contained in the object and its
signed distance from the boundary (positive means external, negative internal).

Only incident intersections are registered, i.e. intersections where the ray
is not tangent to the surface of the object.

Classes
-------
    - ``Plane``:
        A class of hittable objects delimited by a plane.
    - ``Sphere``:
        A class of hittable objects delimited by a sphere.
    - ``HitRecord``:
        A class for recording information about an intersection between a
        hittable and a ray. It includes the point of intersection, the surface
        normal at that point, and the time of the intersection.
"""

from math import sqrt
from math_tools import Vect3D, Ray, dot


class Plane:
    """A class of hittable objects delimited by a plane.

    The normal vector should always be of unit norm. Instance attributes are
    made protected in order to prevent reassignment to a non-unit vector.

    Attributes
    ----------
    _normal : Vect3D
        Orthogonal unit vector, determines the orientation of the plane.
    _shift : float
        Distance of the plane from the one through the origin parallel to it.
    """
    def __init__(self, normal: Vect3D, shift: Vect3D | float):
        """Initializer of the ``Plane`` class.

        An instance is initialized to:
            - The plane of equation ``normal * x == shift`` if ``shift`` is a ``float``; or
            - The unique plane orthogonal to ``normal`` through the point ``shift`` if ``shift`` is a ``Vect3D``.
        Everything is rescaled to ensure the normal vector has unit norm.

        :param Vect3D normal: Normal vector to the surface.
        :param Vect3D | float shift:
        """

        n = normal.norm()
        shift /= n
        normal /= n
        if isinstance(shift, Vect3D):
            shift = dot(normal, shift)
        self._normal = normal
        self._shift = shift

    # Getters
    @property
    def normal(self) -> Vect3D:
        return self._normal

    @property
    def shift(self) -> float:
        return self._shift

    # Methods
    def intercept_time(
            self,
            ray: Ray,
            min_time: float | None = 0,
            max_time: float | None = None,
            tolerance: float = 1.0e-6
    ) -> float | None:
        """Detect time of intersection between ``ray`` and ``self``, if one occurs within an admissible time interval.

        If the line lies entirely on the plane then the intersection is tangential, and it is ignored.

        :param Ray ray: The ray in question.
        :param float | None min_time: Minimum admissible time. ``None`` means no lower bound.
        :param float | None max_time: Maximum admissible time. ``None`` means no upper bound.
        :param float tolerance: Margin of error when estimating incidence between normal and ray direction.
        :return: Time of intersection if one occurs within the specified time interval, ``None`` otherwise.
        """

        t = None
        den = dot(self._normal, ray.direction)
        if abs(den) > tolerance:
            t = (self._shift - dot(self._normal, ray.origin)) / den
        return t if is_within_range(t, min_time, max_time) else None

    def is_hit(
            self,
            ray: Ray,
            min_time: float | None = 0,
            max_time: float | None = None,
            tolerance: float = 1.0e-6
    ) -> bool:
        """Detect whether an intersection between ``self`` and ``ray`` occurs within an admissible time interval.

        :param Ray ray: The ray in question.
        :param float | None min_time: Minimum admissible time. ``None`` means no lower bound.
        :param float | None max_time: Maximum admissible time. ``None`` means no upper bound.
        :param float tolerance: Margin of error when estimating incidence between normal and ray direction.
        :return: ``True`` if an intersection occurs within the specified time interval, ``False`` otherwise.
        """

        return self.intercept_time(ray, min_time, max_time, tolerance) is not None

    def get_hit_record(
            self,
            ray: Ray,
            min_time: float | None = 0,
            max_time: float | None = None,
            tolerance: float = 1.0e-6,
    ) -> 'HitRecord':
        """Record information about the intersection between ``self`` and ``ray``.

        The normal

        :param Ray ray: The ray in question.
        :param float | None min_time: Minimum admissible time. ``None`` means no lower bound.
        :param float | None max_time: Maximum admissible time. ``None`` means no upper bound.
        :param float tolerance: Margin of error when estimating incidence between normal and ray direction.
        :return: A ``HitRecord`` object recording the intersection information.
        """

        t = self.intercept_time(ray, min_time, max_time, tolerance)
        point = ray.point_at(t)
        normal = self._normal if t is not None else None
        return HitRecord(point, normal, t)

    def contains_point(self, point: Vect3D) -> bool:
        """Detect whether ``point`` lies inside the (negative) half-space cut by ``self``."""

        return dot(self._normal, point) >= self._shift

    def distance_from_point(self, point: Vect3D) -> float:
        """Measure signed distance between ``self`` and ``point``."""

        return dot(self._normal, point) - self._shift


class Sphere:
    """A class of hittable objects delimited by a sphere.

    A sphere is characterized by a center and radius, which *may be negative*.
    The sign of the radius indicates whether the solid region delimited
    by the surface is the one inside (positive) or outside (negative) the
    sphere.

    Attributes
    ----------
    _center : Vect3D
    _radius : float
        Positive if the volume of the object is inside the sphere, negative otherwise.
    """

    def __init__(self, center: Vect3D, radius: float):
        """Initializer of the ``Sphere`` class.

        :param Vect3D center: Center of the sphere.
        :param float radius: Radius of the sphere, positive if the volume is inside the sphere, negative otherwise.
        """
        self._center = center
        self._radius = radius

    # Getters
    # @property
    # def center(self) -> Vect3D:
    #     """Getter of the centre (read only)."""
    #
    #     return self._center
    #
    # @property
    # def radius(self) -> float:
    #     """Getter of the radius (read only)."""
    #
    #     return self._radius

    def intercept_time(
            self,
            ray: Ray,
            min_time: float | None = 0,
            max_time: float | None = None,
            tolerance: float = 1.0e-6,
    ) -> float | None:
        """Detect time of intersection between ``ray`` and ``self``, if one occurs within an admissible time interval.

        The intersection condition at time ``t`` reads
        ``(ray.point_at(t) - self._center).norm_squared() == self._radius**2``.
        Since ``ray.point_at(t)`` is linear in ``t``, the above boils down to a quadratic equation
        ``a * t**2 + b * t + c == 0``. A straightforward check shows that the coefficients ``a``,
        ``b``, and ``c`` are as below.

        Existence of real solutions is determined by the discriminant test. If two solutions exist and
        both are within the admissible range, the smaller one (which is reached first by the ray)
        is returned. The case of exactly one real solution, which corresponds to vanishing of the
        discriminant, occurs when the ray is tangent to the sphere, which is not considered a hit.

        :param Ray ray: The ray in question
        :param float | None min_time: Minimum admissible time. ``None`` means no lower bound.
        :param float | None max_time: Maximum admissible time. ``None`` means no upper bound.
        :param float tolerance: Margin of error when estimating the vanishing of the discriminant.
        :return: The smallest of the intersection times of within the specified range, or ``None`` if there are none.
        """

        t = None
        a = ray.direction.norm_squared()
        b_halves = - dot(self._center - ray.origin, ray.direction)
        c = (self._center - ray.origin).norm_squared() - self._radius ** 2
        discriminant = b_halves**2 - a*c
        if discriminant > tolerance / a:
            t = - (b_halves + sqrt(discriminant))/a
            if not is_within_range(t, min_time, max_time):
                t = -2*b_halves - t
            if not is_within_range(t, min_time, max_time):
                t = None
        return t

    def is_hit(
            self,
            ray: Ray,
            min_time: float | None = 0,
            max_time: float | None = None,
            tolerance: float = 1.0e-6,
    ) -> bool:
        """Detect whether an intersection between ``self`` and ``ray`` occurs within an admissible time interval.

        :param Ray ray: The ray in question.
        :param float | None min_time: Minimum admissible time. ``None`` means no lower bound.
        :param float | None max_time: Maximum admissible time. ``None`` means no upper bound.
        :param int tolerance: Margin of error for the discriminant test for tangency.
        :return: ``True`` if an intersection occurs within the specified time interval, ``False`` otherwise.
        """
        
        return self.intercept_time(ray, min_time, max_time, tolerance) is not None

    def get_hit_record(
            self,
            ray: Ray,
            min_time: float | None = 0,
            max_time: float | None = None,
            tolerance: float = 1.0e-6,
    ) -> 'HitRecord':
        """Record information about the intersection between ``self`` and ``ray``.

        :param Ray ray: The ray in question.
        :param float | None min_time: Minimum admissible time. ``None`` means no lower bound.
        :param float | None max_time: Maximum admissible time. ``None`` means no upper bound.
        :param float tolerance: Margin of error for the discriminant test for tangency.
        :return: A ``HitRecord`` object recording the intersection information.
        """
        
        t = self.intercept_time(ray, min_time, max_time, tolerance)
        hit_point = ray.point_at(t)
        normal = (hit_point - self._center) / self._radius if t is not None else None
        return HitRecord(hit_point, normal, t)

    def contains_point(self, point: Vect3D) -> bool:
        """Detect whether ``point`` lies inside the region delimited by the sphere.

        If ``self._radius`` is negative, the region is *outside* the sphere.
        """

        return (point - self._center).norm_squared() / self._radius <= self._radius

    def distance_from_point(self, point: Vect3D):
        """Measure signed distance between ``self`` and ``point``."""

        return (1 if self._radius >= 0 else -1) * (point - self._center).norm() - self._radius


class HitRecord:
    """A class to collect information about the last hit of a ray against the scene.

    Attributes
    ----------
    - point : Vect3D | None
        The point of impact.
    - normal : Vect3D | None
        The surface normal at the point of impact.
    - time : float | None
        The time of impact relative to the parameterization of the ray.
    """

    def __init__(self, point: Vect3D | None, normal: Vect3D | None, time: float | None):
        """Initializer of the ``HitRecord`` class."""

        self.point = point
        self.normal = normal
        self.time = time

    def __str__(self):
        return f"{str(self.point)}, {str(self.normal)}, {str(self.time)}"


def is_within_range(t, min_time=None, max_time=None) -> bool:
    return t is not None and (min_time is None or t >= min_time) and (max_time is None or t <= max_time)
