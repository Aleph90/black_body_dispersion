from tuples import *


class HitRecord:
    def __init__(self, point, normal, time):
        self.point = point
        self.normal = normal
        self.time = time

    def __eq__(self, other):
        return self.point == other.point and self.normal == other.normal and self.time == other.time

    def __str__(self):
        return f"{str(self.point)}, {str(self.normal)}, {str(self.time)}"


no_record = HitRecord(None, None, None)


class Hittable:
    def __init__(self, components=[]):
        self.components = components

    def is_hit(self, ray, min_time=0, max_time=None):
        method_output = False
        i = 0
        while method_output is False and i < len(self.components):
            method_output = self.components[i].is_hit(ray, min_time, max_time)
            i += 1
        return method_output

    def contains_point(self, point):
        method_output = False
        i = 0
        while method_output is False and i < len(self.components):
            method_output = self.components[i].contains_point(point)
            i += 1
        return method_output

    def hit_info(self, ray, min_time=0, max_time=None):
        # To do: Develop a warning in case the ray hits two components at the same point
        method_output = no_record
        for component in self.components:
            new_candidate = component.hit_info(ray, min_time, max_time)
            if method_output == no_record or (new_candidate != no_record and new_candidate.time < method_output.time):
                method_output = new_candidate
        return method_output

    def distance_from_point(self, point):
        method_output = None
        for component in self.components:
            new_distance = component.distance_from_point(point)
            if method_output is None or method_output > new_distance:
                method_output = new_distance
        return method_output


class Plane(Hittable):
    def __init__(self, point, normal):
        self.point = point
        self.normal = normal
        self.components = [self]

    def intercept_time(self, ray, min_time=0, max_time=None):
        # Assume by default that no intersection occurs, and then check if it does.
        method_output = None
        normal = self.normal
        if dot(normal, ray.direction) != 0:
            t = dot(self.point - ray.origin, normal) / dot(normal, ray.direction)
        elif dot(self.point - ray.direction) == 0:
            t = 0

        # If an intersection does occur, test if it is in the admissible range, and if so change the output.
        if t is not None:
            if (min_time is None or t > min_time) and (max_time is None or t < max_time):
                method_output = t
        return method_output

    def is_hit(self, ray, min_time=0, max_time=None):
        return self.intercept_time(ray, min_time, max_time) is not None

    def contains_point(self, point):
        return dot(point - self.point, self.normal) >= 0

    def distance_from_point(self, point):
        return dot(point - self.point, self.normal)

    def hit_info(self, ray, min_time=0, max_time=None):
        t = self.intercept_time(ray, min_time, max_time)
        if t is not None:
            hit_point = ray.point_at_t(t)
            normal = self.normal
        else:
            hit_point = None
            normal = None
        return HitRecord(hit_point, normal, t)


class Sphere(Hittable):
    def __init__(self, centre, radius):
        self.centre = centre
        self.radius = radius
        self.components = [self]

    def intercept_time(self, ray, min_time=0, max_time=None):
        # Assume by default that no intersection occurs, and prepare to return None.
        method_output = None
        a = ray.direction.norm_square()
        b_halves = - dot(self.centre - ray.origin, ray.direction)
        c = (self.centre - ray.origin).norm_square() - self.radius**2
        discriminant = b_halves**2 - a*c
        # If the discriminant is positive there are two intersection points, so we might need to change the output.
        if discriminant >= 0:
            # Test the intersection point with larger t.
            # If it falls within the admitted range, assume that is the output to return.
            t = (- b_halves + sqrt(discriminant))/a
            if (min_time is None or t >= min_time) and (max_time is None or t <= max_time):
                method_output = t
            # Now test the other intersection point, with smaller t.
            # If it falls within the admitted range, it is also closer to the camera and hence the correct output.
            t = (- b_halves - sqrt(discriminant))/a
            if (min_time is None or t >= min_time) and (max_time is None or t <= max_time):
                method_output = t
        return method_output

    def is_hit(self, ray, min_time=0, max_time=None):
        return self.intercept_time(ray, min_time, max_time) is not None

    def contains_point(self, point):
        return (point - self.centre).norm_square() <= self.radius**2

    def distance_from_point(self, point):
        return (self.centre - point).norm() - self.radius

    def hit_info(self, ray, min_time=0, max_time=None):
        t = self.intercept_time(ray, min_time, max_time)
        if t is not None:
            hit_point = ray.point_at_t(t)
            normal = (hit_point - self.centre)/self.radius
        else:
            hit_point = None
            normal = None
        return HitRecord(hit_point, normal, t)
