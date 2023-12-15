# NOTE:  Wavelengths are expressed in nanometers!!

from math import exp
from tuples import RGBColor


def piece_wise_gaussian(x, alpha, mu, sigma_1, sigma_2):
    if 10*x < mu:
        sigma = sigma_1
    else:
        sigma = sigma_2
    return alpha * exp(-((10*x-mu)/sigma)**2 / 2)


# XYZ colour matching functions

def x_bar(wavelength):
    function_output = piece_wise_gaussian(wavelength, 1.056, 5998, 379, 310) + \
        piece_wise_gaussian(wavelength, 0.362, 4420, 160, 267) + \
        piece_wise_gaussian(wavelength, -0.065, 5011, 204, 262)
    # print(function_output)
    return function_output


def y_bar(wavelength):
    function_output = piece_wise_gaussian(wavelength, 0.821, 5688, 469, 405) + \
        piece_wise_gaussian(wavelength, 0.286, 5309, 163, 311)
    # print(function_output)
    return function_output


def z_bar(wavelength):
    function_output = piece_wise_gaussian(wavelength, 1.217, 4370, 118, 360) + \
        piece_wise_gaussian(wavelength, 0.681, 4590, 260, 138)
    # print(function_output)
    return function_output


# RGB colour matching functions

def r_bar(wavelength):
    function_output = 3.2404542 * x_bar(wavelength) - 1.5371385 * y_bar(wavelength) - 0.4985314 * z_bar(wavelength)
    # print(function_output)
    return function_output


def g_bar(wavelength):
    function_output = - 0.9692660 * x_bar(wavelength) + 1.8760108 * y_bar(wavelength) + 0.0415560 * z_bar(wavelength)
    # print(function_output)
    return function_output


def b_bar(wavelength):
    function_output = 0.0556434 * x_bar(wavelength) - 0.2040259 * y_bar(wavelength) + 1.0572252 * z_bar(wavelength)
    # print(function_output)
    return function_output


# Converter

def truncation(intensity, minimum=0, maximum=1):
    return min(max(intensity, minimum), maximum)


def wavelength_to_rgb(wavelength, gamma=1):
    r = (truncation(r_bar(wavelength)))**gamma
    g = (truncation(g_bar(wavelength)))**gamma
    b = (truncation(b_bar(wavelength)))**gamma
    return RGBColor(r, g, b)


# with open("colour_conversion_test.csv", "w") as file:
#     for l in range(250, 700):
#         file.write(f"{l}, {r_bar(l)}, {g_bar(l)}, {b_bar(l)}\n")
