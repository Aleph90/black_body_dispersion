from tuples import *
from random import uniform


# The default value of the positioning depends on the focal length.
# I set the default value to None so the code can detect whether or not an argument is passed by the user.
# If not, a default will be generated later in the method based on the information of the focal length.
#
# If u_hat is given, v_hat will be computed accordingly, and the corresponding argument ignored.
class Camera:
    def __init__(self,
                 image_height,
                 screen_height,
                 aspect_ratio=16.0 / 9.0,
                 focal_length=1,
                 positioning=None,
                 u_hat=None,
                 v_hat=None,
                 anti_aliasing_samples = 100):
        self.aspect_ratio = aspect_ratio
        self.image_height = image_height
        self.image_width = int((image_height * aspect_ratio) // 1)
        self.screen_height = screen_height
        self.screen_width = screen_height * aspect_ratio
        self.pixel_size = self.screen_width/self.image_width

        self.focal_length = focal_length
        if positioning is not None:
            positioning.direction.do_normalise()
            self.positioning = positioning
        else:
            focal_point = Vect(0, -focal_length, 0)
            camera_orientation = Vect(0, 1, 0)
            self.positioning = Ray(focal_point, camera_orientation)

        if u_hat is not None:
            if dot(u_hat, self.positioning.direction) != 0:
                print("The specified screen orientation is not orthogonal to the camera direction!!")
            self.u_hat = u_hat.normalised()
            self.v_hat = cross(self.u_hat, self.positioning.direction)
        elif v_hat is not None:
            if dot(v_hat, self.positioning.direction) != 0:
                print("The specified screen orientation is not orthogonal to the camera direction!!")
            self.v_hat = v_hat.normalised()
            self.u_hat = cross(self.positioning.direction, self.v_hat)
        else:
            # If we are here, it means the screen orientation was not specified.
            # If the camera direction is orthogonal to (1, 0, 0), then that will be the horizontal direction.
            # Otherwise, the default will be the projection of the camera direction to the xz-plane, rotated by -90º.
            if self.positioning.direction.x == 0:
                self.u_hat = Vect(1, 0, 0)
            else:
                self.u_hat = Vect(-self.positioning.direction.z, 0, self.positioning.direction.x).normalised()
            self.v_hat = cross(self.u_hat, self.positioning.direction)

        self.screen_bottom_left = self.positioning.point_at_t(focal_length) \
            - self.u_hat * self.screen_width / 2 - self.v_hat * self.screen_height / 2

        self.anti_aliasing_samples = anti_aliasing_samples

    def focal_point(self):
        return self.positioning.origin

    def random_displacement(self):
        u_displacement = uniform(0, 1)
        while u_displacement == 1:
            u_displacement = uniform(0, 1)
        v_displacement = uniform(0, 1)
        while v_displacement == 1:
            v_displacement = uniform(0, 1)
        return (self.u_hat * u_displacement + self.v_hat * v_displacement) * self.pixel_size
