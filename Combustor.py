"""
Contains the conmbustor heat releasee and heat loss mechanisms
"""
import numpy as np
class Combustor:
    def __init__(self, settings):
        self.Dia = 0.215  # chamber diameter(m)
        self.th_q = 0.0035  # chamber wall thickness(m)
        self.settings = settings
        self.T_exhaust = 0.0
        """self.Re  # nozzle reynolds number
        self.P_rad # kw
        self.P_exh # kW"""
        self.mdot_gas = self.settings.mdot_total # kg / s
        self.cp_gas = 1.3  # kJ / kg - K
        self.sigma = 5.67e-8  # W / m2 - K4
        """self.Nu_gas
        self.Nu_air
        T_quartz_in
        T_quartz_out
        Q_cond
        T_air_main"""


    def heat_loss(self):
        Q_therm = self.settings.P_therm
        Q_in = self.settings.Q_in
        A_comb = 3.14 * self.Dia * 0.49 # m^2
        vol_comb = (3.14/4.0) * (self.Dia **2.0) * 0.49 # m^3
        A_x = 3.14 * (self.Dia ** 2.0) / 4.0
        epsi = 0.05 # 0.165 gas emissivity:0.05
        e_Rad = self.sigma * A_comb * epsi / 1000.0
        e_exh = (self.mdot_gas * self.cp_gas)
        e_rem = (self.mdot_gas * self.cp_gas * self.settings.T_air) + Q_therm + self.sigma * A_comb * epsi * 300.0 ** 4.0 / 1000.0
        coeff = [e_Rad, 0, 0, e_exh, -e_rem]
        roots = np.roots(coeff)
        root_eff = np.max(roots.real[(abs(roots.imag) < 1e-05) & (roots.real > 0.0)])
        self.T_exhaust = root_eff
        print self.T_exhaust
        print self.settings.mdot_total
        print vol_comb
