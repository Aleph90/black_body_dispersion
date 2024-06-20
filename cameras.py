"""This module implements simple cameras for the purpose of ray tracing.

In ray tracing, the color value of each pixel in the output image is computed
by locating the physical position of that pixel inside the camera within
3-dimensional space and retrace the origin of the beam by following its
trajectory in reverse. In a physical camera, the light collected by the
sensor passes through a narrow opening, the objective, which for simplicity
is modelled here as a dimensionless point. The pixels on the receptor, on the
other hand, form a rectangular grid, where each individual pixel has a small
but non-zero extension in space. When back-tracing a beam one needs to decide
from which point exactly to start within the surface area occupied by that
pixel. A common approach is to average over a collection of points drawn
randomly from across the entire area, which also helps improve continuity of
the resulting image.
"""

from math_tools import Vect3D, Matrix3by3, Ray
from random import uniform


class Camera:
    """A class to collect and handle all geometric information of the camera.

    This class stores the geometric information of the camera, including
        - Its location and orientation;
        - Its focal distance, i.e. distance between the film and the objective;
        - Its resolution;
        - The size of the output image.

    All attributes of this class are protected, but getters are provided (read-only).

    The class implements methods for randomly drawing a ray given the indices of a
    pixel on the grid.

    Attributes
    ----------
    _px_size : float
        Physical size of a pixel within the camera.
    _image_width : int
        Width of the output picture, in pixels.
    _image_height : int
        Height of the output picture, in pixels.
    _u_hat : Vect3D
        Frame vector pointing right.
    _v_hat = : Vect3D
        Frame vector pointing down.
    _top_left : Vect3D
        Location of the top left corner of the film.
    _objective : Vect3D
        Location of the objective
    """

    def __init__(
            self,
            image_width: int,
            image_height: int,
            focal_distance: float,
            film_size: float,
            film_dir: str = 'u',
            position: Vect3D = Vect3D(0, 0, 0),
            orientation: Matrix3by3 = Matrix3by3(
                Vect3D(1, 0, 0),
                Vect3D(0, 0, -1),
                Vect3D(0, 1, 0)
            )
    ):
        """Initializer of the Camera class.

        The parameters `image_height` and `image_width`, expressing the size
        of the output image in pixels, are both required. Since the aspect
        ratio of the film inside the camera should match that of the output
        image, only one of its lengths should be specified through
        `film_size`. This is taken by default as width, but it can be changed
        to height using `film_dir`, which can accept 'y', 'z', 'v', 'vertical',
        or 'height', ignoring case.

        The camera is positioned at the point `position` using the middle
        point of the film as anchor. The film is placed on the plane through
        position spanned by the first two columns of `orientation`. The latter
        are also rescaled to the size of a pixel and used as frame vectors
        for the pixel grid, pointing right and down, respectively.
        The objective is aligned with the camera position and placed along
        the *negative* direction of the third column of `orientation` at the
        distance specified by `focal_distance`.

        By default, the camera is placed so the film lies on the `xz`-plane,
        with its sides aligned with the axes and centered at the origin, while
        the objective is placed on the negative `y`-axis.

        **NOTE:** This initializer does *not* check that the given frame is
        non-degenerate, orthonormal, or positively oriented. This allows the
        user to experiment with distortion effects by using a skew frame.
        If the first two columns of `orientation` are linearly dependent,
        then the film degenerates to a slit or a point, and the output image
        becomes a collection of monochromatic lines. If the third column is
        in the span of the first two, then the objective will end up on the
        film and an error will occur when trying to shoot a ray.

        :param int image_width: Width of the output image, in pixels.
        :param int image_height: Height of the output image, in pixels.
        :param float focal_distance: Distance between the film and the objective.
        :param float film_size: Physical size of the film inside the camera.
        :param str film_dir: Specifies whether `film_size`
        :param position:
        :param orientation:
        """

        px_size = film_size / (image_height if film_dir.lower in ['y', 'z', 'v', 'vertical', 'height'] else image_width)
        self._px_size = px_size
        self._image_width = image_width
        self._image_height = image_height
        self._u_hat = orientation[0] * px_size
        self._v_hat = orientation[1] * px_size
        self._top_left = position - image_width * self._u_hat / 2.0 - image_height * self._v_hat / 2.0
        self._objective = position - focal_distance * orientation[2]

    # Getters
    @property
    def px_size(self) -> float:
        """Getter of `_px_size` (read only)."""

        return self._px_size

    @property
    def image_width(self) -> int:
        """Getter of `_image_width` (read only)."""

        return self._image_width

    @property
    def image_height(self) -> int:
        """Getter of `_image_height` (read only)."""

        return self._image_height

    @property
    def u_hat(self) -> Vect3D:
        """Getter of `_u_hat` (read only)."""

        return self._u_hat

    @property
    def v_hat(self) -> Vect3D:
        """Getter of `_v_hat` (read only)."""

        return self._v_hat

    @property
    def top_left(self) -> Vect3D:
        """Getter of `_top_left` (read only)."""

        return self._top_left

    @property
    def objective(self) -> Vect3D:
        """Getter of `_objective` (read only)."""

        return self._objective

    def __getitem__(self, pos: tuple[int, int]) -> Vect3D:
        """Getter operator by a pair of indices.

        Obtain the absolute location of the pixel in position `pos`, by
        its middle point. Indices may be out of range. The position is
        read according to matrix convention (row, column).

        :param int pos: A pair of indices, row than column.
        :return: Absolute position of middle point of pixel in position 'pos'.
        """

        i, j = pos
        return self._top_left + (i+0.5) * self._v_hat + (j+0.5) * self._u_hat

    def get_ray(self, i: int, j: int) -> Ray:
        """Shoot a ray through a random point from the pixel in position (i, j).

        Indices may be out of range. The position is read according to matrix
        convention (row, column).

        :param int i: Row of the pixel.
        :param int j: Column of the pixel.
        :return: Ray shooting from `self._objective` through a random point within the pixel in position (i, j).
        """

        u = j + uniform(0, 1)
        v = j + uniform(0, 1)
        return Ray(self._objective, self._top_left + u*self._u_hat + v*self._v_hat - self._objective)


from math_tools import *
from random import uniform


# The default value of the positioning depends on the focal length.
# The default value is set to None so the code can detect whether an argument is passed by the user.
# If not, a default will be generated later in the method based on the information of the focal length.
#
# If u_hat is given, v_hat will be computed accordingly, and the corresponding argument ignored.
class OldCamera:
    def __init__(self,
                 image_height,
                 film_height,
                 aspect_ratio=16.0 / 9.0,
                 focal_length=1,
                 positioning=None,
                 u_hat=None,
                 v_hat=None,
                 anti_aliasing_samples=100):
        self.aspect_ratio = aspect_ratio
        self.image_height = image_height
        self.image_width = int((image_height * aspect_ratio) // 1)
        self.film_height = film_height
        self.film_width = film_height * aspect_ratio
        self.pixel_size = self.film_width / self.image_width
        self.focal_length = focal_length

        if positioning is not None:
            positioning.direction.do_normalise()
            self.positioning = positioning
        else:
            focal_point = Vect3D(0, -focal_length, 0)
            camera_orientation = Vect3D(0, 1, 0)
            self.positioning = Ray(focal_point, camera_orientation)

        if u_hat is not None:
            if dot(u_hat, self.positioning.direction) != 0:
                print("The specified screen orientation is not orthogonal to the camera direction!!")
            self.u_hat = u_hat.normalized()
            self.v_hat = cross(self.u_hat, self.positioning.direction)
        elif v_hat is not None:
            if dot(v_hat, self.positioning.direction) != 0:
                print("The specified screen orientation is not orthogonal to the camera direction!!")
            self.v_hat = v_hat.normalized()
            self.u_hat = cross(self.positioning.direction, self.v_hat)
        else:
            # If we are here, it means the screen orientation was not specified.
            # If the camera direction is orthogonal to (1, 0, 0), then that will be the horizontal direction.
            # Otherwise, the default will be the projection of the camera direction to the xz-plane, rotated by -90º.
            if self.positioning.direction.x == 0:
                self.u_hat = Vect3D(1, 0, 0)
            else:
                self.u_hat = Vect3D(-self.positioning.direction.z, 0, self.positioning.direction.x).normalized()
            self.v_hat = cross(self.u_hat, self.positioning.direction)

        self.screen_bottom_left = self.positioning.point_at_t(focal_length) \
                                  - self.u_hat * self.film_width / 2 - self.v_hat * self.film_height / 2

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
