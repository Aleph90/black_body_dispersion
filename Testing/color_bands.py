from color_utilities import *

image = 'P3\n180\n390\n255\n'
for freq in range(400, 790):
    wl = 3e5/freq
    color = wl2ciergb(wl).round_components(scale=255)
    image += f'{color.r} {color.g} {color.b}\n'*30
    color = wl2ciergb(wl).apply_gamma(0.4).round_components(scale=255)
    image += f'{color.r} {color.g} {color.b}\n'*30
    color = wl2lsrgb(wl).round_components(scale=255)
    image += f'{color.r} {color.g} {color.b}\n'*30
    color = wl2srgb(wl).round_components(scale=255)
    image += f'{color.r} {color.g} {color.b}\n'*30
    color = wl2argb(wl).round_components(scale=255)
    image += f'{color.r} {color.g} {color.b}\n'*30
    color = wl2argb(wl).apply_gamma(0.4).round_components(scale=255)
    image += f'{color.r} {color.g} {color.b}\n'*30

with open('color_bands.ppm', 'w') as output_file:
    output_file.write(image)
