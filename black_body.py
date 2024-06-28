"""A module to implement Planck's law for the radiance of a black body.

Planck's law states that, under appropriate conditions, a body at a set
temperature `T` (typically expressed in Kelvin) emits energy in form
of electromagnetic radiation, and provides an explicit formula for its
spectral radiance, i.e. the amount of power emitted per unit area,
solid angle, and wavelength/frequency. This formula is a function of
a wavelength `w` and the temperature of the body, and it reads

`I(w, T) = N w^{-5} / (e^{A T/w} - 1)`

where the parameters `N` (for "normalization") and `B` are combinations
of physical constants such as the Boltzmann and Planck constant and the
speed of light.

There are at least two useful interpretations of this function when
studying the light emitted by the body along a certain ray:
    - As the intensity of the light in each wavelength;
    - As the probability distribution that a beam has a given wavelength.
Correspondingly, this module develops two classes built on these
respective interpretations.
"""

from math import exp, pi
from random import choice

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
    _wl_scale: float
        Scaling factor to account for different choices of unit for wavelength
        out of meters, micrometers, and nanometers.
    _exp_coeff: float
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
        self._wl_scale = label2scale[units]
        self._exp_coeff = PLANCK * SPEED_OF_LIGHT / (BOLTZMANN * temperature * self._wl_scale)
        self._amplitude = 2 * PLANCK * SPEED_OF_LIGHT / (self._wl_scale ** 5)

    # Getters
    @property
    def temperature(self) -> float:
        return self._temperature

    @property
    def wl_scale(self) -> str:
        return scale2label[self._wl_scale]

    @property
    def exp_coefficient(self) -> str:
        return self._exp_coeff

    @property
    def amplitude(self) -> float:
        return self._amplitude

    def __call__(self, wl: float) -> float:
        """Call method of the ``BlackBodyIntensity`` class.

        Takes a wavelength as argument, expressed in meters, micrometers, or
        nanometers in accordance with the value of the ``_wl_scale``
        attribute, and computes the intensity of the emitted radiation in that
        wavelength, *up to rescaling* by ``_amplitude``.

        :param float wl: Input wavelength.
        :return: Intensity of emitted radiation in the wavelength ``wl``.
        """

        return 1 / ((exp(self._exp_coeff / wl) - 1) * wl ** 5)


class BlackBodyDistribution:
    """Probability distribution of wavelengths given the temperature.

    Sample random wavelengths from a discrete approximation of the Planck
    distribution. The values are drawn from a predefined range and spaced out by
    intervals of equal total probability. In other words, the approximate
    distribution is constructed by mapping a uniform probability on a discrete
    set onto the Planck distribution, restricted to the specified interval.
    The theoretical exact map that achieves this in the continuum is the quantile
    function, i.e. the inverse of the antiderivative of the density function.
    This, however, is not readily available as a closed formula, so the class
    is based on numerical approximations and a tabular representation instead.

    An initial cutoff for the variable rho may be provided as an optional
    argument, as values of rho below 1/800 cause erros in the computation of the
    density function, which however takes negligible values in that range (under
    1.0e-33 at rho = 1/100.

    Many of the attributes of this class are used only upon instantiation and
    stored for debugging purposes.

    Attributes
    ----------
    _table : list[int]
        Table of sample wavelengths, spaced at intervals of even probability.
        Choosing an index uniformly gives an approximation of drawing from
        the Planck distribution.
    _temperature : float
        Temperature of the body modelled by the instance.
    _wl_scale : float
        Indicates the scale of measurement of wavelengths (meters, micrometers,
        nanometers). It is expressed as the ratio between the unit in use and
        meters (m -> 1.0, mum -> 1.0e-6, nm -> 1.0e-9).
    _low_wl : float
        The lower end of the allowed range of wavelengths.
    _high_wl : float
        The upper end of the allowed range of wavelengths.
    _steps : int
        The number of intervals into which the range is divided. This is one
        *less* than the total number of samples.
    _delta : float
        The width of the rectangles used in the integral approximation.
    _start_rho : float
        Starting point for the integration variable rho.
    """

    def __init__(
            self,
            temperature: float,
            interval: tuple[float, float],
            steps: int,
            delta: float,
            units: str | float = 'nm',
            start_rho: float = 0.005,
    ):
        """Initializer of the BlackBodyDistribution class.

        Generate and store a list of wavelengths from a prescribed range, located at
        intervals of even probability according to the Planck law. The inputs specify
        physical properties of the body (temperature), the shape of the distribution,
        and some controls for technical aspects of the execution. They are mostly never
        used again after initialization but are stored as read-only attributes for
        debugging purposes.

        The procedure uses numerical approximation of integrals via the rectangles method.
        Two different changes of variables, rho and tau, are used in the respective regimes
        of small and large wavelengths. This improves efficiency by reducing both ends
        of the integration domain to bounded intervals, and by simplifying the expression
        of the integrand to one that can be computed using fewer operations. Also in the
        interest of efficiency, the calculation uses a version of the density that is *not*
        normalized; however, its mass can be computed analytically, and a single one-time
        rescaling on the output allows saving a much larger number of operations.

        :param float temperature: Temperature of the black body.
        :param tuple[float, float] interval: Pair of floats delimiting the range of wavelengths to draw from.
        :param int steps: Number of sub-intervals.
        :param float delta: Width of the rectangle used for integral approximation.
        :param str | float units: Specifies the unit of measurement of wavelengths (meters, micrometers, nanometers).
        """

        self._temperature = temperature
        self._low_wl = interval[0]
        self._high_wl = interval[1]
        self._steps = steps
        self._delta = delta
        self._wl_scale = label2scale[units]
        self._start_rho = start_rho

        # Lookup table
        self._table = []

        # Utility parameters and functions
        a = PLANCK * SPEED_OF_LIGHT / (BOLTZMANN * temperature * self._wl_scale)

        def density_by_rho(x: float):
            return 1 / ((exp(1/x) - 1) * x**5)

        def density_by_tau(x: float):
            return (x**3) / (exp(x) - 1)

        # Determine the total probability "mass" of the lower tail. Start from 0.0 and
        # then add mass in approximate increments using the variable rho = wl/a. Stop
        # when wl reaches the lower end of the interval.
        u = 0.0
        rho = start_rho
        max_rho = self._low_wl / a
        while rho < max_rho:
            rho += delta
            u += delta * density_by_rho(rho)

        # Cut out the total probability "mass" of the upper tail. The total mass of this
        # un-normalized distribution is pi**4 / 15, so start from that value and subtract
        # mass in approximate increments using the variable tau = a/wl, so the interval
        # becomes bounded between 0.0 and a/self._high_wl. This quantity will later be used
        # as a stop condition, so it is called cap.
        cap = pi ** 4 / 15.0
        tau = 0.0
        max_tau = a / self._high_wl
        while tau < max_tau:
            tau += delta
            cap -= delta * density_by_tau(tau)

        # The remaining mass needs to be subdivided into chunks of equal weight. Starting from
        # the lower end of the wavelength interval, use the variable tau again to add increments
        # of mass to the variable u. Every time a chunk of appropriate size is formed, record the
        # current wavelength in the lookup table.
        step = (cap - u) / steps
        flag = u
        tau = 1/rho
        while len(self._table) <= self._steps:
            if u >= flag:
                self._table.append(a/tau)
                flag += step
            tau -= delta
            u += delta * density_by_tau(tau)

    # Getters
    @property
    def table(self) -> list[float]:
        return self._table

    @property
    def temperature(self) -> float:
        return self._temperature

    @property
    def scale(self) -> float:
        return self._wl_scale

    @property
    def low_wl(self) -> float:
        return self._low_wl

    @property
    def high_wl(self) -> float:
        return self._high_wl

    @property
    def steps(self) -> int:
        return self._steps

    @property
    def delta(self) -> float:
        return self._delta

    # Draw
    def draw(self) -> float:
        return choice(self._table)
