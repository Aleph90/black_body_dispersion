"""A module for materials and related classes.

Materials in this model are responsible for determining how a surface
interacts with light --- whether it emits, absorbs, reflects, or refracts
it. This module implements several materials as classes, each including
a method ``scatter`` which determines if and how light is deflected after
hitting a surface, given the essential information about the geometry of
the impact.

The parameters of some of these phenomena can sometimes depend on the
wavelength of the light involved. Such is the case for the refractive index
of some materials, which can be approximated by the Sellmeier equation. It
takes the form
.. math::
    \\eta^{2} = 1 + \\sum \\frac{B_{i} \\lambda^{2}}{\\lambda^{2} - C_{i}}

where :math:`\\lambda` is the wavelength and the :math:`B_{i}`'s and
:math:`C_{i}`'s are appropriate coefficients. This module includes classes
which implement single components of this formula and the full
wavelength-to-index function.
"""

from math import sqrt
from math_tools import Ray, SumFun
from beams import *
from hittables import HitRecord
from black_body import label2scale


class SellmeierComponent:
    """A class for individual components of the Sellmeier equation.

    All attributes are protected. The coefficients should not be modified
    manually. A change in the scale should be reflected in the parameter ``c``,
    so the scale should also not be changed manually but only through the
    provided method.

    Attributes
    ----------
    _b : float
        Coefficient at the numerator.
    _c : float
        Coefficient at the denominator.
    _scale: float
        Specifies the unit of length for wavelengths.
    """

    def __init__(self, b: float, c: float, units: float | str = 'nm'):
        """Initializer of the ``SellmeierComponent`` class.

        Stores the three coefficients appearing in the component, as well as
        a specification of the unit of measurement. The coefficient `c` is
        adjusted to the scale.

        :param float b: Coefficient at the numerator.
        :param float c: Coefficient at the denominator.
        :param float | str units: Specifies the unit of measurement for wavelengths.
        """

        self._scale = label2scale[units]
        self._b = b
        self._c = c / (self._scale ** 2)

    # Getters
    @property
    def scale(self) -> float:
        return self._scale

    @property
    def b(self) -> float:
        return self._b

    @property
    def c(self) -> float:
        return self._c

    def rescale(self, new_units: float | str) -> None:
        """Change the scale of lengths used for wavelengths.

        :param float | str new_units: The new unit.
        :return: ``None``
        """
        self._c *= self._scale**2
        self._scale = label2scale[new_units]
        self._c /= self.scale**2

    def __call__(self, wl: float) -> float:
        """Call method.

        :param float wl: A wavelength, assumed expressed in the same units as ``self._scale``.
        :return: The contribution given by the component at the specified wavelength.
        """
        return self.b / (1 - self.c / wl**2)


class Sellmeier:
    """A class for the full Sellmeier equation of a material.

    The Sellmeier equation consists of several components. This class collects them
    all together into a single object.
    """
    def __init__(self, *args):
        comps = []
        units = set()
        for arg in args:
            comps.append(arg if isinstance(arg, SellmeierComponent) else SellmeierComponent(arg))
            units.add(comps[-1].scale)
        self._scale = units.pop() if len(units) == 1 else label2scale['nm']
        if units:
            for comp in comps:
                comp.rescale(self._scale)
        self._correction = SumFun(*comps)

    @property
    def scale(self) -> float:
        return self._scale

    @property
    def correction(self) -> SumFun:
        return self._correction

    def __call__(self, wl: float) -> float:
        return sqrt(1 + self._correction(wl))


class Material:
    def __init__(self):
        pass

    def scatter(self, beam: Beam, rec: HitRecord):
        raise NotImplementedError('Material not implemented.')


class Dielectric(Material):
    """A class to model materials featuring refraction.

    Attributes
    ----------
    _refractive_idx : float
        The refractive index of the material inside the surface (read only).
    """

    def __init__(self, refractive_idx: float | None = None):
        """Initializer of the ``Dielectric`` class.

        :param refractive_idx: Refractive index of the material.
        """

        super().__init__()
        self._refractive_idx = refractive_idx

    @property
    def refractive_idx(self) -> float:
        return self._refractive_idx

    def get_refractive_idx(self, beam: Beam) -> float:
        """Compute the refractive index.

        The ray argument, while redundant for this implementation,
        is included to maintain consistency with subclasses.

        :param Beam beam: Not used, see above.
        :return: The refractive index of the material.
        """

        return self._refractive_idx

    def scatter(self, beam: Beam, rec: HitRecord):
        if beam.direction.dot_with(rec.normal) < 0:
            # The beam is propagating against the surface normal, i.e.
            # from outside the hittable object into it. The refractive
            # index before refraction is the last in the beam's record,
            # after refraction it is the material's. Since the beam is
            # entering the surface, the new index needs to be recorded.
            eta_before = beam.ref_rec.last()
            eta_after = beam.ref_rec.append(self.get_refractive_idx(beam))
        else:
            # The beam is propagating along the surface normal, i.e.
            # from inside it toward the outside. The refractive index
            # before refraction is the material's, which is also the
            # last index recorded by the beam. Since beam is leaving the
            # material, this index needs to be removed from the record.
            # The index after refraction is the second last, which after
            # popping the "before" index is simply the last.
            eta_before = beam.ref_rec.pop()
            eta_after = beam.ref_rec.last()
        beam.direction = beam.direction.refract(rec.normal, eta_before, eta_after)
        beam.origin = rec.point
        return True


class SpectralDielectric(Dielectric):
    """A class of dielectric materials with refractive index dependent of wavelength.

    Uses the same method as the parent class to scatter a beam, but overrides the
    ``get_refractive_idx`` method to fetch the refractive index from the wavelength
    of the beam.

    Attributes
    ----------
    _sellmeier : Sellmeier
        Function that provides the refractive index in terms of wavelength.
    """

    def __init__(self, sellmeier: Sellmeier, *args, **kwargs):
        """Initializer of the ``SpectralDielectric`` class.

        Requires the data of a ``Sellmeier`` object.

        :param Sellmeier sellmeier: Sellmeier data.
        :param args: Arguments for the parent initializer.
        :param kwargs: Keyword arguments for the parent initializer.
        """

        super().__init__(*args, **kwargs)
        self._sellmeier = sellmeier

    def get_refractive_idx(self, beam: SpectralBeam) -> float:
        """Compute refractive index.

        The conversion from wavelength to refractive index is carried
        out by ``self._sellmeier``, while the information of the
        wavelength is stored in the beam. Requires a *spectral* beam.

        :param SpectralBeam beam: A spectral beam carrying the wavelength information.
        :return: The refractive index of the material for the wavelength of the specified beam.
        """

        return self._sellmeier(beam.wl)


class LightSource(Material):
    """A class of materials that emit light

    """
