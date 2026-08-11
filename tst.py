from sample_mass_calcs.xas_sample import XRaySample

tst = XRaySample("Cu", "Cu", "K", mu_total=1,
                 density=8.96,area = 1.3)

tst.calculate_thickness()
tst.length_unit = "um"
tst.mass_unit = "mg"

print(f"THICKNESS = {tst.thickness}, ABS STEP = {tst.absorbance_step}")

tst.calculate_mass()
