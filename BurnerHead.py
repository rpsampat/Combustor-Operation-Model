import numpy as np
import cantera as ct


class BurnerHead:
    def __init__(self, settings, T_heater):
        self.N = 12.0  # number of nozzles
        self.nozzle_dia = 0.008  # nozzle diameter(m)
        self.nozzle_dia_fuel = 0.001  # fuel nozzle diameter(m)
        self.st = settings
        self.vel_nozzle = 0  # velocity at nozzle (m/s)
        self.vel_nozzle_fuel = 0  # velocity at fuel nozzle (m/s)
        self.Re = 0  # nozzle Reynolds number
        self.Re_fuel = 0  # fuel nozzle Reynolds number
        self.T_air = T_heater
        self.flowrates()

    def flowrates(self):
        area_nozzle = (np.pi/4.0)*self.nozzle_dia**2.0
        area_nozzle_fuel = (np.pi/4.0)*self.nozzle_dia_fuel**2.0
        air = ct.Solution('air.cti')
        air.TP = self.T_air, ct.one_atm
        gas = ct.Solution('gri30.cti')
        gas.TPX = 293.15, ct.one_atm, {'CH4': 1.0}
        self.vel_nozzle = self.st.mdot_nozzle/(air.density*area_nozzle*self.N)
        self.vel_nozzle_fuel = self.st.mdot_nozzle_fuel / (gas.density * area_nozzle_fuel * self.N)
        self.Re = air.density*self.vel_nozzle*self.nozzle_dia/air.viscosity
        self.Re_fuel = gas.density * self.vel_nozzle_fuel * self.nozzle_dia_fuel / gas.viscosity
