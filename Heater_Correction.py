import os
import math
import Thermocouple as TC

class Heater_Correction:
    def __init__(self):
        self.dia_elem = 0.025  # heater element diameter
        self.dia_pipe = 0.044  # heater total pipe diameter
        self.mdot_air_total = 1000.0  # lnpm
        self.elem_num = 3.0  # number of heater elements

    def thermocouple_corr_control(self):
        tc = TC.Thermocouple()
        tc.dia_tc = 0.001  # thermocouple diameter
        tc.mdot_cool = self.mdot_air_total/self.elem_num  # air flow rate (lnpm)
        tc.T_air = 600 + 273.15  # K
        T_wall = 100 + 273.15 # K
        tc.area_conv = (math.pi/4.0) * self.dia_elem ** 2.0
        root_eff, Q_rad, Q_rad_chamber = tc.heat_balance(T_wall,area_cond=0.0)
        print root_eff-273.15

    def thermocouple_corr_measure(self):
        tc = TC.Thermocouple()
        tc.dia_tc = 0.002  # thermocouple diameter
        tc.mdot_cool = self.mdot_air_total  # air flow rate (lnpm)
        tc.T_air = 600 + 273.15  # K
        T_wall = 100 + 273.15 # K
        tc.area_conv = (math.pi/4.0) * self.dia_pipe ** 2.0
        root_eff, Q_rad, Q_rad_chamber = tc.heat_balance(T_wall, area_cond=0.0)
        print root_eff-273.15

if __name__ == "__main__":
    hc = Heater_Correction()
    hc.thermocouple_corr_control()
    hc.thermocouple_corr_measure()
