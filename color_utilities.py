"""A module to handle color and conversion operations.

The main purpose of this module is to implement RGB color  and matching
functions, i.e. the conversion from wavelength to RGB representation of
monochromatic light.

Color representation is implemented as a subclass ``Color`` of ``Vect3D``.
It overrides some of the parent methods to ensure consistency of return
types and provides aliases ``r``, ``g``, and ``b`` for the components. It
also implements component-by-component transformations, such as gamma
correction and clipping of values within a prescribed range, and export
methods to save color as pixels in a .ppm file.

Color conversion for monochromatic light from wavelength to RGB is
implemented by using a different color space, XYZ, as an intermediate
step. The matching functions for this space, i.e. the converters from
pure wavelength to XYZ representation, can be approximated explicitly
as sums of piece-wise Gaussian functions, i.e. functions defined
independently on two separate intervals by Gaussian expressions with
different parameters to the left and right of the peak.

There are several versions of RGB. This module implements conversion
functions from the XYZ color space to
    - CIE 1931 RGB, obtained from XYZ by a linear transformation;
    - Linear sRGB, also linear in XYZ;
    - sRGB, or standard RGB, developed by HP and Microsoft for use
        with screens and printers. It is obtained from linear sRGB
        by a non-linear transformation on the three components.
    - Adobe RGB, which offers a broader gamut of colors compared to
        sRGB. It is also linear in XYZ, possibly up to a gamma
        correction of coefficient 256/263 ~ 0.9734.


Classes
-------
    - ``Color``:
        RGB representation of color, subclass of ``Vect3D``. Overrides some
        of the ``Vect3D`` and implements new ones for component-wise
        transformations, plus some utilities for saving to .ppm format.
    - ``PWGaussian``:
        Piece-wise Gaussian function.

Functions
---------
    - ``wl2xyz``:
        Full XYZ color matching function. Converts wavelength to XYZ.

    - ``wl2ciergb``:
        CIE RGB color matching function. Converts wavelength to CIE RGB.

    - ``wl2lsrgb``:
        Linear standard RGB color matching function. Converts wavelength
        to linear standard RGB.

    - ``wl2srgb``:
        Standard RGB color matching function. Converts wavelength to
        standard RGB.

    - ``wl2argb``:
        Adobe RGB color matching function. Converts wavelength to Adobe RGB.

Other objects
-------------
    - ``PWGaussian wl2[xyz][012]``:
        Individual summands for the components of the XYZ color matching function.
    - ``SumFun wl2[xyz]``:
        Components of the XYZ color matching function.
    - ``Matrix3by3 xyz2ciergb1``:
        Matrix for the XYZ to CIE RGB transformation.
    - ``Matrix3by3 xyz2srgb``:
        Matrix for the XYZ to standard RGB transformation.
    - ``Matrix3by3 xyz2argb``:
        Matrix for the XYZ to Adobe RGB transformation.

NOTE: All wavelengths are expressed in *nanometers*.
"""

from math import exp
from math_tools import Vect3D, Matrix3by3, SumFun


# Color
class Color(Vect3D):
    """Custom class to represent color in RGB format.

    These objects are linear in nature, and many of the operations involving them
    are similar to those of 3-dimensional vectors, including algebraic operations
    like sum, rescaling, linear combinations, and matrix multiplication. As such,
    the class is implemented as a subclass of ``Vect3D``. Many of the method
    definitions are overridden to ensure that the objects returned by all
    operations are of the correct type. In addition, property methods are defined
    to give aliases ``r`` (red), ``g`` (green), and ``b`` (blue) to the ``x``-,
    ``y``-, and ``z``-components of a color expressed as a vector.
    """

    def __init__(self, r: float, g: float, b: float):
        """Initializer of the ``RGBColor`` class.

        :param float r: Red component.
        :param float g: Green component.
        :param float b: Blue component.
        """

        super().__init__(r, g, b)

    def __eq__(self, other) -> bool:
        """Equality operator for the ``RGBColor`` class.

        Overrides the parent method by adding the condition that both
        terms be instances of the ``RGBColor`` class.

        :param Color other: The color to be compared with ``self``.
        :return: ``True`` if both ``self`` and ``other`` are instances of ``RGBColor``
            and are equal as ``Vect3D`` objects.
        """

        return isinstance(other, Color) and super().__eq__(other)

    # Representations and type conversion
    def __repr__(self) -> str:
        """Override representation operator to reflect the new type.

        :return: A string representation of ``self``, specifying the ``RGBColor`` type.
        """

        return f'RGBColor{tuple(self._entries)}'

    def print_ppm_line(self) -> str:
        """Generate a line that can be added to a ppm file to represent a pixel.

        :return: A string representing the color as a pixel in ppm format.
        """

        return f'{self.r} {self.g} {self.b}\n'

    def round_components(self, scale: int = 1):
        """Produce a rescaled version of ``self`` with integer components.

        :param int scale: Rescaling factor, default to 1, i.e. no rescaling.
        :return: A new RGB color with values rescaled by ``scale`` and integer components.
        """

        return Color(*[int((scale * e) // 1.0) for e in self._entries])

    # Access
    @property
    def r(self) -> float:
        """Red component, alias for ``x``.

        :return: The value of the red component of a color.
        """

        return self.x

    @property
    def g(self) -> float:
        """Green component, alias for ``y``.

        :return: The value of the green component of a color.
        """

        return self.y

    @property
    def b(self) -> float:
        """Blue component, alias for ``z``.

        :return: The value of the blue component of a color.
        """

        return self.z

    @r.setter
    def r(self, value: float) -> None:
        """Set the value of the red component.

        :param float value: Value to assign to the red component.
        :return: ``None``.
        """

        self.x = value

    @g.setter
    def g(self, value: float) -> None:
        """Set the value of the green component.

        :param float value: Value to assign to the green component.
        :return: ``None``.
        """

        self.y = value

    @b.setter
    def b(self, value: float) -> None:
        """Set the value of the blue component.

        :param float value: Value to assign to the blue component.
        :return: ``None``.
        """

        self.z = value

    # Algebra
    def __add__(self, other: 'Color') -> 'Color':
        """Addition operator for the ``RGBColor`` class.

        Overridden to ensure that the returned object is of the correct type, i.e.
        `RGBColor` if both summands are instances of this class. Since ``self``
        certainly is, the return type should be the same as that of the other summand.

        :param Color other: The color to be added to ``self``.
        :return: A new color representing the superposition of ``self`` and ``other``.
        """

        return self.cast(super().__add__(other))

    def __sub__(self, other: 'Color') -> 'Color':
        """Subtraction operator for the ``RGBColor`` class.

        Overridden to ensure that the returned object is of the correct type, i.e.
        ``RGBColor`` if both members are instances of this class. Since ``self``
        certainly is, the return type should be the same as that of the other term.

        :param Color other: The color to be subtracted from ``self``.
        :return: A new color representing the difference of ``self`` and ``other``.
        """

        return self.cast(super().__sub__(other))

    def __rmul__(self, scalar: float) -> 'Color':
        """Rescaling operator for the ``RGBColor`` class.

        Overridden to ensure that the returned object is of the correct type, i.e. ``RGBColor``.

        :param float scalar: The factor by which to rescale ``self``.
        :return: A new color representing the rescaled color.
        """

        return self.cast(super().__rmul__(scalar))

    def __mul__(self, scalar: float) -> 'Color':
        """Rescaling operator for the ``RGBColor`` class by multiplication from the right.

        Overridden to ensure that the returned object is of the correct type, i.e. ``RGBColor``.

        :param float scalar: The factor by which to rescale ``self``.
        :return: A new color representing the rescaled color.
        """

        return self.cast(super().__mul__(scalar))

    def __truediv__(self, scalar: float) -> 'Color':
        """Rescaling operator for the ``RGBColor`` class.

        Overridden to ensure that the returned object is of the correct type, i.e. ``RGBColor``.

        :param float scalar: The factor by which to rescale ``self``.
        :return: A new color representing the rescaled color.
        """

        return self.cast(super().__truediv__(scalar))

    def __neg__(self) -> 'Color':
        """Negation operator for the ``RGBColor`` class.

        May not be meaningful for colors, but it is convenient to define subtraction.
        Overridden to ensure the correct return type, i.e. ``RGBColor``.

        :return: A new "color" with components opposite to those of ``self``.
        """
        return self.cast(super().__neg__())

    # Color transformations.
    def transform_components(self, f) -> 'Color':
        """Apply a transformation to all components individually.

        :param f: Transformation to apply to all components.
        :return: A new RGB color with all components transformed according to ``f``.
        """

        return Color(*[f(e) for e in self._entries])

    def clip(self, minimum: float = 0, maximum: float = 1):
        """Clipping function to ensure color is within allowable range.

        If any component is outside the allowable range it will be replaced
        with the value of closer end.

        :param minimum: Lower end of allowable range.
        :param maximum: Upper end of allowable range.
        :return: A new ``RGBColor`` object with all three components clipped within range.
        """

        return self.transform_components(lambda x: min(max(x, minimum), maximum))

    def apply_gamma(self, gamma: float):
        """Gamma correction.

        Gamma correction amounts to raising all components to a certain power,
        the gamma coefficient.
        When gamma < 1, the effect of this operation is to increase all values,
        resulting in lighter colors. It also enhances the separation between
        darker colors and attenuates it for lighter ones.
        When gamma > 1, the effect is the opposite: all values are decreased and
        colors made darker, with separation enhanced in the lighter range.

        The three RGB components are assumed to be expressed as floating
        point numbers between 0.0 and 1.0.

        :param float gamma: Coefficient for gamma correction.
        :return: A new RGB color with gamma correction applied.
        """

        return self.transform_components(lambda x: x**gamma)


class PWGaussian:
    """Piece-wise Gaussian function.

    This function is defined by two different formulas, or lobes,
    for values of ``x`` before and after a certain threshold ``mu``.
    Both lobes are Gaussian functions with the same maximum value
    `alpha` attained at ``mu``, and only differ by their "width".
    """

    def __init__(self, alpha: float, mu: float, tau1: float, tau2: float):
        """Initializer of the ``PWGaussian`` class.

        :param float alpha: Amplitude of the function, i.e. its maximum value.
        :param float mu: Location of the peak.
        :param float tau1: Inverse spread of left lobe.
        :param float tau2: Inverse spread of the right lobe.
        """
        self._alpha = alpha
        self._mu = mu
        self._tau1 = tau1
        self._tau2 = tau2

    def __call__(self, wl: float) -> float:
        tau = self._tau1 if wl < self._mu else self._tau2
        return self._alpha * exp(-(tau * (wl - self._mu)) ** 2 / 2)


# Piece-wise Gaussian functions used to approximate color matching functions.
wl2x0 = PWGaussian(0.362, 442.0, 0.0624, 0.0374)
wl2x1 = PWGaussian(1.056, 599.8, 0.0264, 0.0323)
wl2x2 = PWGaussian(-0.065, 501.1, 0.0490, 0.0382)
wl2y0 = PWGaussian(0.821, 568.8, 0.0213, 0.0247)
wl2y1 = PWGaussian(0.286, 530.9, 0.0613, 0.0322)
wl2z0 = PWGaussian(1.217, 437.0, 0.0845, 0.0278)
wl2z1 = PWGaussian(0.681, 459.0, 0.0385, 0.0725)

# XYZ Color matching functions as single components.
wl2x = SumFun(wl2x0, wl2x1, wl2x2)
wl2y = SumFun(wl2y0, wl2y1)
wl2z = SumFun(wl2z0, wl2z1)


def wl2xyz(wl: float) -> Vect3D:
    """XYZ color matching function.

    :param float wl: Wavelength of the monochromatic light to be converted.
    :return: XYZ representation of monochromatic light of wavelength ``wl`` as a ``Vect3D`` object.
    """

    return Vect3D(wl2x(wl), wl2y(wl), wl2z(wl))


# Matrix of the XYZ to RGB linear map.
xyz2ciergb = Matrix3by3(
    Color(+2.36461385, -0.51516621, +0.0052037),
    Color(-0.89654057, +1.4264081, -0.01440816),
    Color(-0.46807328, +0.0887581, +1.00920446)
)

xyz2srgb = Matrix3by3(
    Color(+3.2406, -0.9689, +0.0557),
    Color(-1.5372, +1.8758, -0.2040),
    Color(-0.4986, +0.0415, +1.0570)
)

xyz2argb = Matrix3by3(
    Color(+2.04159, -0.96924, +0.01344),
    Color(-0.56501, +1.87598, -0.11836),
    Color(-0.34473, +0.04156, +1.01517)
)


def wl2ciergb(wl: float) -> Color:
    """CIE RGB color matching function.

    :param float wl: Wavelength of monochromatic light.
    :return: RGB representation of monochromatic light of wavelength ``wl`` according to the CIE standard.
    """
    return (xyz2ciergb * wl2xyz(wl)).clip()


def wl2lsrgb(wl: float) -> Color:
    """Linear sRGB color matching function.

    :param float wl: wavelength of monochromatic light.
    :return: RGB representation of monochromatic light of wavelength ``wl`` according to the sRGB standard,
        without gamma correction.
    """

    return (xyz2srgb * wl2xyz(wl)).clip()


def wl2srgb(wl: float) -> Color:
    """sRGB color matching function.

    :param float wl: wavelength of monochromatic light.
    :return: RGB representation of monochromatic light of wavelength ``wl`` according to the sRGB standard.
    """

    return wl2lsrgb(wl).transform_components(lambda x: 12.92*x if x <= 0.0031308 else 1.055*x**0.4166 - 0.055)


def wl2argb(wl: float) -> Color:
    """Adobe RGB color matching function.

    :param float wl: wavelength of monochromatic light.
    :return: RGB representation of monochromatic light of wavelength ``wl`` according to the Adobe RGB standard.
    """

    return (xyz2argb * wl2xyz(wl)).clip()
