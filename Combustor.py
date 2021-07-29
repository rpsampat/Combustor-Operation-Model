class Combustor:
    def __init__(self):
        self.Dia = 0.215 # chamber diameter(m)
        self.th_q = 0.0035 # chamber wall thickness(m)
        self.v_gas
        self.T_exhaust
        self.T_ad
        self.Re #nozzle reynolds number
        self.P_rad # kw
        self.P_exh # kW
        self.mdot_gas # kg / s
        self.cp_gas = 1.3; # kJ / kg - K
        self.sigma = 5.67e-8; # W / m2 - K4
        self.Nu_gas
        self.Nu_air
        T_quartz_in
        T_quartz_out
        Q_cond
        T_air_main

    def heat_loss(self):
        