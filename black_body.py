"""A module to implement Planck's law for the radiance of a black body.

Planck's law states that, under appropriate conditions, a body at a set
temperature `T` (typically expressed in Kelvin) emits energy in form
of electromagnetic radiation, and provides an explicit formula for its
spectral radiance, i.e. the amount of power emitted per unit area,
solid angle, and wavelength/frequency. This formula is a function of
a wavelength `w` and the temperature of the body, and it reads

`I(w, T) = N w^{-5} / (e^{A T/w} - 1)`

where the parameters `N` (for "normalization" and `B` are combinations
of physical constants such as the Boltzmann and Planck constant and the
speed of light.

There are at least two useful interpretations of this function when
studying the light emitted by the body along a certain ray:
    - As the intensity of the light in each wavelength;
    - As the probability distribution that a beam has a given wavelength.
Correspondingly, this module develops two classes built on these
respective interpretations.
"""

from math import exp

# Physical constants
BOLTZMANN = 1.380649e-23        # m^2 kg / s^2 K
PLANCK = 6.62607015e-34         # m^2 kg / s
SPEED_OF_LIGHT = 2.99792458e8   # m / s

# Scale conversion
"""Wavelength may be expressed in meters, micrometers, or nanometers.
This dict will be of use later to switch from one scale to another.
"""
label2scale = {
    'meters': 1.0,
    'metres': 1.0,
    'meter': 1.0,
    'metre': 1.0,
    'm': 1.0,
    1.0: 1.0,
    0: 1.0,

    'micrometers': 1.0e-6,
    'micrometres': 1.0e-6,
    'micrometer': 1.0e-6,
    'micrometre': 1.0e-6,
    'micro': 1.0e-6,
    'mum': 1.0e-6,
    'mu': 1.0e-6,
    1.0e-6: 1.0e-6,
    -6: 1.0e-6,

    'nanometers': 1.0e-9,
    'nanometres': 1.0e-9,
    'nanometer': 1.0e-9,
    'nanometre': 1.0e-9,
    'nano': 1.0e-9,
    'nm': 1.0e-9,
    'n': 1.0e-9,
    1.0e-9: 1.0e-9,
    -9: 1.0e-9,
}

scale2label = {
    1.0: 'm',
    1.0e-6: 'mum',
    1.0e-9: 'nm',
}


class BlackBodyIntensity:
    """Radiance of a black body of given temperature as a function of wavelength.

    The class provides a ``__call__`` method that computes the value of the
    radiance at a given wavelength. Because the function is to be called a large
    number of times throughout execution, all attributes are read-only as a way
    to avoid inconsistencies. As a further way to improve efficiency, the constant
    coefficients involved in the function are computed upon initialization and
    stored as attributes.
    The function involves a constant overall multiplicative factor. In practice,
    however, the function outputs are processed linearly and eventually rescaled
    to adjust the picture brightness, so applying the multiplicative factor on
    function call is inefficient and its effect eventually overwritten anyway.
    For that reason, the constant factor is stored as an attribute but not applied
    by the call method.

    Attributes
    ----------
    _temperature: float
        Temperature of the black body, expressed in Kelvin.
    _wavelength_scale: float
        Scaling factor to account for different choices of unit for wavelength
        out of meters, micrometers, and nanometers.
    _exp_coefficient: float
        Multiplicative coefficient inside the exponential term of the radiance
        function. Stored as an attribute to save operations.
    _amplitude: float
        Overall multiplicative factor, stored but not applied upon call.
    """

    def __init__(self, temperature: float, units: str | float = 'nm'):
        """Initializer of the ``BlackBodyIntensity`` class.

        Stores the temperature of the black body as well as all the necessary
        coefficients as read-only attributes. Also records a scale parameter
        to account for the different possible choices of unit of measurement
        of the input wavelength. The accepted inputs and corresponding scales
        are stored in the module attribute ``label2scale``, and include.

        :param float temperature: Temperature of the black body.
        :param units: Specifies the units of the input wavelengths.
        """

        self._temperature = temperature
        self._wavelength_scale = label2scale[units]
        self._exp_coefficient = PLANCK * SPEED_OF_LIGHT / (BOLTZMANN * temperature * self._wavelength_scale)
        self._amplitude = 2 * PLANCK * SPEED_OF_LIGHT / (self._wavelength_scale ** 5)

    # Getters
    @property
    def temperature(self) -> float:
        return self._temperature

    @property
    def wavelength_scale(self) -> str:
        return scale2label[self._wavelength_scale]

    @property
    def exp_coefficient(self) -> str:
        return self._exp_coefficient

    @property
    def amplitude(self) -> float:
        return self._amplitude

    def __call__(self, wl: float) -> float:
        """Call method of the ``BlackBodyIntensity`` class.

        Takes a wavelength as argument, expressed in meters, micrometers, or
        nanometers in accordance with the value of the ``_wavelength_scale``
        attribute, and computes the intensity of the emitted radiation in that
        wavelength, *up to rescaling* by ``_amplitude``.

        :param float wl: Input wavelength.
        :return: Intensity of emitted radiation in the wavelength ``wl``.
        """

        return 1 / ((exp(self._exp_coefficient/wl) - 1) * wl**5)
