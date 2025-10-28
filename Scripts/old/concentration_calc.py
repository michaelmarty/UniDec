import numpy as np
import scipy.constants as const

na=const.physical_constants["Avogadro constant"][0]
pi=const.pi

def sphere(radius):
    volume=4/3. * pi *(radius)**3.
    area=4 * pi * (radius)**2.
    return volume, area

def cube(width):
    volume = width**3.
    area = 6 * width**2.
    return volume, area

def minimum_molar_box(conc):
    nm3=pow(10.,24.)/(na*conc)
    nm=pow(nm3,1/3.)
    return nm

print(minimum_molar_box(0.1))
print([minimum_molar_box(i) for i in [1,0.1,0.01,0.001,0.0001,0.00001,0.000001,0.0000001]])