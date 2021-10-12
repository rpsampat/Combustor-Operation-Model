import os
import math
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits import axisartist
class Thermocouple:
    def __init__(self):
        self.epsi_tc = 0.79  # thermocouple emissivity
        self.epsi_wall = 0.79  # combustor wall emissivity
        self.sigma = 5.6704e-8  # stefan-boltzman constant
        self.mdot_cool = 3000  # wall cooling flow (lnpm)
        self.T_air = 60+273.15  # K
        self.dia_tc = 0.001  # thermocouple diameter (m)
        self.dia_chamber = 0.215  # combustion chamber diameter (m)
        self.Pr_air = 0.7
        self.k_steel = 15  # W/mK
        self.k_tc_wall = self.k_steel  # conductivity of thermocouple - wall interface
        self.area_conv = 0  # thermocouple convection effective area

    def conv_area_calc_ext(self):
        """
        Effective thermocouple convection area for thermocouples located outside the combustion chamber
        perpendicular to the wall
        :return:
        """
        area = (math.pi / 4.0) * ((0.215 + 0.09) ** 2 - (0.215 ** 2))  # m^2

        return area

    def convection_thermocouple(self):
        """
        Convection between thermocouple and external cooling flow modelled as cross flow around cylinder
        :return:
        """
        air = ct.Solution('air.cti')
        air.TP = self.T_air, ct.one_atm
        mdot = self.mdot_cool *1.225/60000  #kg/s
        vel_tc = mdot/(air.density*self.area_conv)
        Re_tc = air.density*vel_tc*self.dia_tc/air.viscosity
        Nu_lam = 0.664*(Re_tc**0.5)*(self.Pr_air**0.333)
        Nu_turb = 0.037*(Re_tc**0.8)*self.Pr_air/(1+2.443*(Re_tc**-0.1)*(self.Pr_air**(2.0/3.0)-1))
        Nu = 0.3 + math.sqrt(Nu_lam**2+Nu_turb**2)
        h = Nu*air.thermal_conductivity/self.dia_tc
        #print "Nu=", Nu

        return h

    def heat_balance(self, T_wall, area_cond):
        """
        A T^4 + B T + C = 0
        :return:
        """
        dx = 0.001
        h = self.convection_thermocouple()
        area_rad_tc = (math.pi/4.0)*self.dia_tc**2 + math.pi*self.dia_tc*0.10
        area_rad_chamber = math.pi*self.dia_chamber*0.1

        area_conv = math.pi*self.dia_tc*0.10
        A = -1.0* self.epsi_tc *self.sigma*area_rad_tc
        B = (-self.k_steel/dx)*area_cond- h*area_conv
        C = self.epsi_wall*self.sigma*T_wall**4.0*area_rad_tc + (self.k_steel/dx)*T_wall*area_cond \
            + h * self.T_air * area_conv
        coeff = [A, 0.0, 0.0, B, C]
        roots = np.roots(coeff)
        root_eff = roots.real[(abs(roots.imag)< 1e-05) & (roots.real>0.0)]
        Q_rad = -1.0*self.epsi_tc *self.sigma*area_rad_tc*root_eff[0]**4.0\
                + self.epsi_wall*self.sigma*T_wall**4.0*area_rad_tc
        Q_rad_chamber = self.epsi_wall*self.sigma*T_wall**4.0*area_rad_chamber
        return root_eff, Q_rad, Q_rad_chamber

    def main(self, path):
        self.area_conv = self.conv_area_calc_ext()
        area_cond = (math.pi / 4.0) * self.dia_tc ** 2
        T_wall = 200.0 + 100.0*np.array(range(20))  # C
        T_wall_k = T_wall +273.15 # K
        T_eff = np.zeros(len(T_wall_k))
        Q_rad = np.zeros(len(T_wall_k))
        Q_rad_chamb = np.zeros(len(T_wall_k))
        for i in range(len(T_wall_k)):
            roots, rad, rad_chamb = self.heat_balance(T_wall_k[i], area_cond)
            T_eff[i] = roots[0]-273.15
            Q_rad[i] = rad
            Q_rad_chamb[i] = rad_chamb/1000
            #print(T_wall_k[i])

        fig, ax = plt.subplots()
        ax.plot(T_eff, label = "Effective Temp")
        ax.plot(T_wall, label = "Wall Temp")
        ax.set_xlabel("Index")
        ax.set_ylabel("Temperature (C)")
        ax.legend()
        plt.grid(True)
        plt.savefig(path + '/Temperature_comparison', dpi=300, bbox_inches='tight')
        plt.close()

        fig, ax = plt.subplots()
        ax.plot(T_wall, T_eff, label="Effective Temp")
        ax.set_xlabel("Wall Temperture (C)")
        ax.set_ylabel("Effective Temperature (C)")
        ax.legend()
        plt.grid(True)
        plt.savefig(path + '/Effective_Temp', dpi=300, bbox_inches='tight')
        plt.close()

        fig2 = plt.figure()
        host = host_subplot(111, axes_class=axisartist.Axes)
        # plt.subplots_adjust(right=0.75)
        ax1 = host.twinx()
        p, = host.plot(T_wall, Q_rad, label="Radiation on TC")
        p1, = ax1.plot(T_wall, Q_rad_chamb, label="Radiation from Chamber")
        ax1.axis["right"].toggle(all=True)
        host.set_xlabel("Temperature")
        host.set_ylabel("Net Radiation on TC (W)")
        ax1.set_ylabel("Radiation from Chamber (kW)")
        host.legend()
        plt.grid(True)
        plt.savefig(path + '/Heatloss', dpi=300, bbox_inches='tight')
        plt.close()
        plt.show()


if __name__ == "__main__":
    TC = Thermocouple()
    path = os.getcwd()
    TC.main(path)




