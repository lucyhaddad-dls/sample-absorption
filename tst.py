
from src.xas_sample import XRaySample


tst = XRaySample(formula="TiO", absorber="Ti",
                 edge="K", mass=1)


# print(f"initial mass unit={tst.mass.unit}, density={tst.density.unit}")

# # change mass to be in mg.
# tst.mass.unit = "mg"
# print(f"next: mass={tst.mass.unit}, density={tst.density.unit}")


print(tst.mass_unit, tst.length_unit)
print(tst.mass)
tst.mass_unit = "mg"

# print(tst.mass_unit, tst.length_unit)
# print(tst.mass)