import numpy as np
import Settings as st
import BurnerHead as bh
import Combustor as comb



class Run:

    def __init__(self):
        self.P_min = 60
        self.P_max = 70
        self.p_step = 10
        self.P_Range = len(np.arange(self.P_min, self.P_max,self.p_step))

        self.eq_min = 0.8
        self.eq_max = 0.9
        self.eq_step = 0.1
        self.eq_range = len(np.arange(self.eq_min, self.eq_max, self.eq_step))

        self.P_therm = np.zeros((self.eq_range, self.P_Range))
        self.P_rad = np.zeros((self.eq_range, self.P_Range))
        self.P_cond = np.zeros((self.eq_range, self.P_Range))
        self.T_gas = np.zeros((self.eq_range, self.P_Range))
        self.T_gas_ad = np.zeros((self.eq_range, self.P_Range))
        self.T_quartz = np.zeros((self.eq_range, self.P_Range))
        self.vdot_air = np.zeros((self.eq_range, self.P_Range))
        self.vdot_air_cool = np.zeros((self.eq_range, self.P_Range))
        self.v_nozzle = np.zeros((self.eq_range, self.P_Range))
        self.v_nozzle_fuel = np.zeros((self.eq_range, self.P_Range))
        self.volflow_fuel = np.zeros((self.eq_range, self.P_Range))
        self.v_gas = np.zeros((self.eq_range, self.P_Range))
        self.T_heater = np.zeros((self.eq_range, self.P_Range))
        self.T_exhaust = np.zeros((self.eq_range, self.P_Range))
        num_equi = 0;
        O2_perc = 20.95;
        self.fuel = 'CH4'

    def const_dilution(self):
        """
        Combustor run at constant external dilution but varying equivalence ratio and thermal power
        :return:
        """
        phi = self.eq_min
        for i in range(self.eq_range):
            power = self.P_min
            for j in range(self.P_Range):
                print "Power=",power
                print "Phi=",phi
                T_heater = 273.15+400
                settings = st.Settings(self.fuel, power, phi, T_heater)
                burner = bh.BurnerHead(settings, T_heater)
                combustor = comb.Combustor(settings)
                combustor.combustor()

                power = power+self.p_step

            phi = phi + self.eq_step

if __name__ == "__main__":
    run = Run()
    run.const_dilution()