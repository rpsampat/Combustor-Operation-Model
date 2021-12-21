"""
Contains the combustor heat release and heat loss mechanisms
"""
import numpy as np
import Combustor_CRN as ccrn
import CRNPIV_data as crnpiv

class Combustor:
    def __init__(self, settings):
        self.Dia = 0.215  # chamber diameter(m)
        self.th_q = 0.0035  # chamber wall thickness(m)
        self.settings = settings
        self.T_exhaust = 0.0
        self.vol_comb = (3.14/4.0) * (self.Dia **2.0) * 0.49 # m^3
        self.A_comb = 3.14 * self.Dia * 0.49  # m^2
        """self.Re  # nozzle reynolds number
        self.P_rad # kw
        self.P_exh # kW"""
        self.mdot_gas = self.settings.mdot_total # kg / s
        self.cp_gas = 1.3  # kJ / kg - K
        self.sigma = 5.67e-8  # W / m2 - K4
        self.adiabatic = True
        """self.Nu_gas
        self.Nu_air
        T_quartz_in
        T_quartz_out
        Q_cond
        T_air_main"""

    def heat_loss(self):
        """
        Calculates radiative heat loss through a transparent quartz combustor
        :return:
        Q_rad:
        """
        Q_therm = self.settings.P_therm
        Q_in = self.settings.Q_in
        A_x = 3.14 * (self.Dia ** 2.0) / 4.0
        epsi = 0.05 # 0.165 gas emissivity:0.05
        e_Rad = self.sigma * self.A_comb * epsi / 1000.0
        e_exh = (self.mdot_gas * self.cp_gas)
        e_rem = (self.mdot_gas * self.cp_gas * self.settings.T_air) + Q_therm \
                + self.sigma * self.A_comb * epsi * 300.0 ** 4.0 / 1000.0
        coeff = [e_Rad, 0, 0, e_exh, -e_rem]
        roots = np.roots(coeff)
        root_eff = np.max(roots.real[(abs(roots.imag) < 1e-05) & (roots.real > 0.0)])
        self.T_exhaust = root_eff
        Q_rad = e_Rad*self.T_exhaust**4.0

        return Q_rad

    def combustor_chemistry(self, Q_rad,case_num):
        """
        Calls the Combustor_CRN class which contains a CRN model of the combustor. This function also takes heat loss
        as an input and passes it to the CRN.
        :param Q_rad:
        :return:
        """
        U = Q_rad/self.A_comb
        crn = ccrn.CRN(self.settings,U,self.A_comb)
        piv_dat = crnpiv.CRNPIV_data()
        piv_dat.main(case_num)
        crn.connect = piv_dat.connect
        crn.alpha_mat = piv_dat.alpha_mat
        crn.crn_def = piv_dat.reactors
        crn.outlet = piv_dat.outlet
        crn.vol_comb = piv_dat.vol_comb_sector
        crn.wall_reactor = piv_dat.reactor_wall
        crn.adiabatic = self.adiabatic
        crn.case_num = case_num
        if case_num==3:
            crn.heat_loss_reduc_fact =0.80
        else:
            crn.heat_loss_reduc_fact =1.0
        print "Volume comb =",crn.vol_comb
        #t, n, c = crn.main()

        data = crn.main2()

        return data

    def combustor(self,case_num):
        qrad = self.heat_loss()
        #t,n,c = self.combustor_chemistry(qrad)
        data_set = self.combustor_chemistry(0,case_num)

        return data_set



