
from src.xas_sample import XRaySample
import matplotlib.pyplot as plt

tst = XRaySample(formula="TiO", absorber="Ti",
                 edge="K", mass=1, thickness=20,
                 density=1, mass_unit="kg",
                 energy_unit="eV")

element = tst.elements[0]
element.mass_unit = "lb"

fig, ax = plt.subplots()
for element in tst.elements:
    ax.plot(tst.energy.value, element.mass_absorption.value,
            label=element.mass_unit)
ax.legend()
plt.show()