# NOTE:  Wavelengths are expressed in nanometers!!

from math import exp
from tuples import RGBColor


class PieceWiseGaussian:
    def __init__(self, alpha: float, mu: int, sigma_1: int, sigma_2: int):
        self.alpha = alpha
        self.mu = mu
        self.sigma_1 = sigma_1
        self.sigma_2 = sigma_2

    def __call__(self, x: float) -> float:
        sigma = self.sigma_1 if 10*x < self.mu else self.sigma_2
        return self.alpha * exp(-((10*x - self.mu)/sigma)**2 / 2)


class XYZColorMatching:
    def __init__(self, components: list[PieceWiseGaussian]):
        self.components = components

    def __call__(self, x: float) -> float:
        return sum([component(x) for component in self.components])


class XYZtoRGBConverter:
    def __init__(self, c_x: float, c_y: float, c_z: float):
        self.c_x = c_x
        self.c_y = c_y
        self.c_z = c_z

    def __call__(self, x: float, y: float, z: float) -> float:
        return self.c_x * x + self.c_y * y + self.c_z * z


# XYZ color matching function
wavelength_to_x_bar = XYZColorMatching([PieceWiseGaussian(1.056, 5998, 379, 310),
                                        PieceWiseGaussian(0.362, 4420, 160, 267),
                                        PieceWiseGaussian(-0.065, 5011, 204, 262)])
wavelength_to_y_bar = XYZColorMatching([PieceWiseGaussian(0.821, 5688, 469, 405),
                                        PieceWiseGaussian(0.286, 5309, 163, 311)])
wavelength_to_z_bar = XYZColorMatching([PieceWiseGaussian(1.217, 4370, 118, 360),
                                        PieceWiseGaussian(0.681, 4590, 260, 138)])


# RGB color matching functions
r_bar = XYZtoRGBConverter(3.2404542, -1.5371385, -0.4985314)
g_bar = XYZtoRGBConverter(-0.9692660, 1.8760108, 0.0415560)
b_bar = XYZtoRGBConverter(0.0556434, -0.2040259, 1.0572252)


# Converter
def truncation(intensity, minimum=0, maximum=1):
    return min(max(intensity, minimum), maximum)


def wavelength_to_rgb(wl: float, gamma=1) -> RGBColor:
    x = wavelength_to_x_bar(wl)
    y = wavelength_to_x_bar(wl)
    z = wavelength_to_z_bar(wl)
    r = truncation(r_bar(x, y, z))**gamma
    g = truncation(g_bar(x, y, z))**gamma
    b = truncation(b_bar(x, y, z))**gamma
    return RGBColor(r, g, b)
