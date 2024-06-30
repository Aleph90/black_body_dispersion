"""A simple module for enriched versions of the ``Ray`` class.

This module implements subclasses of ``Ray`` that hold information
about a ray as modelling physical light. This includes a record of
the refractive indices of the materials traversed by the ray and,
in the spectral version, the wavelength.

Classes
-------
RefractiveRecord:
    Maintains a record of all refractive indices of the materials
    entered by a beam of light. Implemented as a stack: an index
    is added to the record when the beam enters the material and
    removed when it leaves.
Beam (Ray)
    Ray with the additional information of a ``RefractiveRecord``.
SpectralBeam (Beam)
    Beam with the additional information of a wavelength. Represents
    a beam of monochromatic light.
"""

from math_tools import Vect3D, Ray


class RefractiveRecord:
    """A record of refractive indices of materials entered by a beam.

    In our model each hittable object comes with a material representing
    the properties of its interior. Rather than assigning both interior
    and exterior materials to objects and dealing with the resulting
    potential inconsistencies, we equip beams with the information of
    previously encountered materials. Beams are initialized with the
    record of all materials of the objects enclosing the camera, and the
    record is automatically updated as a stack whenever they cross a
    surface.

    Attributes
    ----------
    _rec : list[float]
        Represents the list of refractive indices encountered up to a
        given point. Can be accessed but not modified manually.
    """
    def __init__(self, rec: list[float] | None = None):
        """Initializer of the ``RefractiveRecord`` class.

        :param rec: List of pre-existing refractive indices, or ``None``.
            If ``None`` is provided, the stack is initialized empty.
        """

        self._rec = rec if rec is not None else []

    @property
    def rec(self) -> list[float]:
        return self._rec

    def last(self) -> float:
        """Access the last entered material without modifying the record.
        If the record is empty, assume the beam is in the vacuum and return
        ``1.0``.

        :return: The refractive index of the last material entered, or `None`.
        """

        return self._rec[-1] if self._rec else 1.0

    def pop(self) -> float:
        """Access the last entered material and remove it from the record.
        If the record is empty, assume the beam is in the vacuum and return
        ``1.0``.

        :return: The refractive index of the last material entered.
        """

        return self._rec.pop() if self._rec else 1.0

    def append(self, idx: float) -> float:
        """Adds the refractive index of a new material to the top of the stack.

        :param float idx: The new index to be added to the stack.
        :return: The added index.
        """

        self._rec.append(idx)
        return idx


class Beam(Ray):
    """A class to model a physical beam of light with the additional
    information of the materials surrounding it.
    """

    def __init__(
            self,
            origin: Vect3D,
            direction: Vect3D,
            ref_rec: RefractiveRecord | None = None,
    ):
        """Initializer of the ``Beam`` class.

        :param Vect3D origin: Starting point of the beam.
        :param Vect3D direction: Direction vector of the beam.
        :param RefractiveRecord ref_rec: A record of the materials
            surrounding the source that emitted the beam.
        """

        super().__init__(origin, direction)
        self.ref_rec = ref_rec if ref_rec is not None else RefractiveRecord()


class SpectralBeam(Beam):
    """A class to model a beam of monochromatic light of given wavelength.
    """

    def __init__(
            self,
            origin: Vect3D,
            direction: Vect3D,
            wavelength: float,
            ref_rec: RefractiveRecord = RefractiveRecord(),
    ):
        """Initializer of the ``SpectralBeam`` class.

        :param Vect3D origin: Starting point of the beam.
        :param Vect3D direction: Direction vector of the beam.
        :param wavelength: Wavelength of the beam.
        :param RefractiveRecord ref_rec: A record of the materials
            surrounding the source that emitted the beam.
        """
        super().__init__(origin, direction, ref_rec)
        self._wl = wavelength

    # Getters
    @property
    def wl(self) -> float:
        return self._wl
