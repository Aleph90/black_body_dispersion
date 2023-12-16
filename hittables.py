from tuples import *


def is_within_range(t, min_time=None, max_time=None) -> bool:
    return t is not None and (min_time is None or t >= min_time) and (max_time is None or t <= max_time)


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


class CompositeHittable:
    def __init__(self, components=None):
        self.components = components if components else []

    def is_hit(self, ray, min_time=0, max_time=None):
        hit_detected = False
        i = 0
        while not hit_detected and i < len(self.components):
            hit_detected = self.components[i].is_hit(ray, min_time, max_time)
            i += 1
        return hit_detected

    def contains_point(self, point):
        container_found = False
        i = 0
        while not container_found and i < len(self.components):
            container_found = self.components[i].contains_point(point)
            i += 1
        return container_found

    def hit_info(self, ray, min_time=0, max_time=None):
        # To do: Develop a warning in case the ray hits two components at the same point
        record = no_record
        for component in self.components:
            new_candidate = component.hit_info(ray, min_time, max_time)
            if record == no_record or (new_candidate != no_record and new_candidate.time < record.time):
                record = new_candidate
        return record

    def distance_from_point(self, point):
        distance = None
        for component in self.components:
            new_distance = component.distance_from_point(point)
            if distance is None or new_distance < distance:
                distance = new_distance
        return distance


class Plane:
    def __init__(self, point, normal):
        self.point = point
        self.normal = normal
        self.components = [self]

    def intercept_time(self, ray, min_time=0, max_time=None):
        # Assume by default that no intersection occurs, and then check if it does.
        t = None
        if dot(self.normal, ray.direction) != 0:
            t = dot(self.point - ray.origin, self.normal) / dot(self.normal, ray.direction)
        elif dot(self.point - ray.origin, self.normal) == 0:
            t = 0
        return t if is_within_range(t, min_time, max_time) else None

    def is_hit(self, ray, min_time=0, max_time=None):
        return self.intercept_time(ray, min_time, max_time) is not None

    def contains_point(self, point):
        return dot(point - self.point, self.normal) >= 0

    def distance_from_point(self, point):
        return dot(point - self.point, self.normal)

    def hit_info(self, ray, min_time=0, max_time=None):
        t = self.intercept_time(ray, min_time, max_time)
        return HitRecord(ray.point_at_t(t), self.normal, t) if t is not None else no_record


class Sphere:
    def __init__(self, center, radius):
        self.center = center
        self.radius = radius
        self.components = [self]

    def intercept_time(self, ray, min_time=0, max_time=None):
        # Assume by default that no intersection occurs, and prepare to return None.
        t = None
        a = ray.direction.norm_square()
        b_halves = - dot(self.center - ray.origin, ray.direction)
        c = (self.center - ray.origin).norm_square() - self.radius ** 2
        discriminant = b_halves**2 - a*c
        # If the discriminant is positive there are two intersection points.
        # We want the earliest admissible solution, or None if neither is.
        # If the smaller solution is admissible, we will accept it.
        # Otherwise, we test the other one.
        if discriminant >= 0:
            t = (- b_halves - sqrt(discriminant))/a
            if not is_within_range(t, min_time, max_time):
                t = (- b_halves + sqrt(discriminant))/a
            if not is_within_range(t, min_time, max_time):
                t = None
        return t

    def is_hit(self, ray, min_time=0, max_time=None):
        return self.intercept_time(ray, min_time, max_time) is not None

    def contains_point(self, point):
        return (point - self.center).norm_square() <= self.radius**2

    def distance_from_point(self, point):
        return (self.center - point).norm() - self.radius

    def hit_info(self, ray, min_time=0, max_time=None):
        t = self.intercept_time(ray, min_time, max_time)
        if t is not None:
            hit_point = ray.point_at_t(t)
            normal = (hit_point - self.center) / self.radius
        else:
            hit_point = None
            normal = None
        return HitRecord(hit_point, normal, t)
