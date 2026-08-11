from dataclasses import dataclass
from pint import Unit, UnitRegistry
import numpy as np

units = UnitRegistry(system='mks')

units.define("eV = [energy] = EV = ev = electronvolt")
units.define("GeV = 1e3 eV = GEV = gev = gigaelectronvolt")


@dataclass
class Measurement:
    """
    A measurement has an attribute `value` and property\
    `unit`. \\
    
    Example:
        ```
        velocity = Measurement(value = [0, 1, 2], 
                            _unit = units.m / units.s**2)

        velocity.to(units.cm / units.s**2)
        ```
    """
    value: int | float | np.ndarray | None
    _unit: Unit | str

    @property
    def unit(self)->None:
        """
        Set unit value.
        """
        return self._unit

    @unit.setter
    def unit(self, value:str|Unit):
        if isinstance(value, Unit):
            self._unit = value
        elif isinstance(value, str):
            value = getattr(units, value)
            self._unit = value

    def to(self, unit:str|Unit):
        """
        Convert value to new unit.
        """
        if isinstance(unit, Unit):
            unit = str(unit)

        if self.value is not None:
            val = self.value * self.unit
            new = val.to(unit)
            self.value = new._magnitude
        self.unit = getattr(units, unit)

    def __mul__(self, other):
        return Measurement(value = self.value * other.value,
                          _unit = self.unit*other.unit)
    
    def __add__(self, other):
        if self.unit != other.unit:
            try:
                other.to(self.unit)
            except:
                raise ValueError(f"Units {self.unit} {other.unit} \
                                 not able to be summed.")

        return Measurement(value = self.value + other.value,
                           _unit = self.unit)

    def __sub__(self, other):
        if self.unit != other.unit:
            try:
                other.to(self.unit)
            except:
                raise ValueError(f"Units {self.unit} {other.unit} \
                                    not able to be subtracted.")
        return Measurement(value = self.value - other.value,
                           _unit = self.unit)

    def __truediv__(self, other):
        if self.unit != other.unit:
            try:
                # this will return a scalar - no longer a measurement.
                other.to(self.unit)
                return self.value/other.value

            except:
                return Measurement(value = self.value/other.value,
                                   _unit = self.unit/other.unit)

    def __pow__(self, exp):
        return Measurement(value= self.value**exp, 
                           _unit = self.unit**exp)


class SampleMeasurements:
    "Work in progress.."
    def __init__(self,
                density:float|int|None = None,
                surface_density:float|int|None = None,
                mass:float|int|None = None,
                area:float|int|None = None,
                volume:float|int|None = None,
                thickness:float|int|None = None,
                mass_unit:Unit|str = "g",
                length_unit:Unit|str = "cm",
                energy_unit:Unit|str = "gev"):
        """
        Test to see if sample measurements should be their own
        object type..
        """

        # density is the important measurement for rescaling mass-absorption coef.

        if isinstance(mass_unit, Unit):
            self._mass_unit = mass_unit
        elif isinstance(mass_unit, str):
            self._mass_unit = getattr(units, mass_unit)

        if isinstance(length_unit, Unit):
            self._length_unit = length_unit
        elif isinstance(length_unit, str):
            self._length_unit = getattr(units, length_unit)


        if isinstance(energy_unit, Unit):
            self._energy_unit = energy_unit
        elif isinstance(energy_unit, str):
            self._energy_unit = getattr(units, energy_unit)


        self._thickness = Measurement(value = thickness, _unit = self.length_unit)
        self._area = Measurement(value = area, _unit = self.length_unit**2)
        self._volume = Measurement(value = volume, _unit = self.length_unit**3)

        self._mass = Measurement(value = mass, _unit=self.mass_unit)
        self._surface_density = Measurement(value = surface_density, _unit = self.mass_unit/self.length_unit**2)
        self._density = Measurement(value = density, _unit = self.mass_unit/self.length_unit**3)
    

    @property
    def mass_unit(self):
        return self._mass_unit

    @mass_unit.setter
    def mass_unit(self, val:Unit|str):
        if isinstance(val, str):
            self._mass_unit = getattr(units, val)
        else:
            self._mass_unit = val

        # convert all mass-dependent values:
    
        return self._mass_unit
    
    @property
    def length_unit(self):
        return self._length_unit

    @length_unit.setter
    def length_unit(self, val:Unit|str):
        if isinstance(val, Unit):
            self._length_unit = val
        else:
            self._length_unit = getattr(units, val)

        # convert all length-dependent values:

        return self._length_unit

    @property
    def energy_unit(self):
        return self._energy_unit

    @energy_unit.setter
    def energy_unit(self, val:str|Unit):
        if isinstance(val, Unit):
            self._energy_unit = val
        else:
            self._energy_unit = getattr(units, val)

        self.edge_energy.to(self.energy_unit)
        self.energy.to(self.energy_unit)

  