import cantera as ct


class Settings:

    def __init__(self, fuel, P_therm, phi, T_heater,vdot_N2):
        self.P_therm = P_therm  # kW
        self.fuel = fuel  # 'CH4'/'H2'/'CH4+H2'
        self.LHV = 0
        self.phi = phi
        self.AF_stoic = 0
        self.O2_perc = 0
        self.mdot_fuel = 0  # kg/s
        self.mdot_air = 0 # kg/s
        self.mdot_N2 = 0 # kg/s
        self.mdot_total = 0 #kg/s
        self.vdot_fuel = 0 # lnpm
        self.vdot_air = 0 # lnpm
        self.vdot_N2 =vdot_N2 # lnpm
        self.T_air = T_heater
        self.T_fuel = 273.15 + 15
        self.T_in = 0
        self.Q_in =0
        self.Y_comp = {}
        self.enthalpy_inlet = 0
        self.enthalpy_outlet =0
        self.main()


    def LHV_CH4_calc(self, gas):
        # CH4 + 2 O2 = CO2 + 2 H2O
        gas.TPX = 298, ct.one_atm, 'CH4:1, O2:2'
        h1 = gas.enthalpy_mass
        Y_CH4 = gas['CH4'].Y[0]  # returns an array, of which we only want the first element

        # set state to complete combustion products without changing T or P
        gas.TPX = None, None, 'CO2:1, H2O:2'
        h2 = gas.enthalpy_mass
        lhv = -(h2 - h1) / Y_CH4 / 1e6
        #print('LHV = {:.3f} MJ/kg'.format(-(h2 - h1) / Y_CH4 / 1e6))

        return lhv

    def LHV_H2_calc(self, gas):
        # 2 H2 + O2 = 2 H2O
        gas.TPX = 298, ct.one_atm, 'H2:2, O2:1'
        h1 = gas.enthalpy_mass
        Y_H2 = gas['H2'].Y[0]  # returns an array, of which we only want the first element

        # set state to complete combustion products without changing T or P
        gas.TPX = None, None, 'H2O:2'
        h2 = gas.enthalpy_mass
        lhv = -(h2 - h1) / Y_H2 / 1e6
        #print('LHV = {:.3f} MJ/kg'.format(-(h2 - h1) / Y_H2 / 1e6))

        return lhv

    def main(self):
        gas = ct.Solution('gri30.cti')
        gas_diluent = ct.Solution('gri30.cti')
        air = ct.Solution('air.cti')
        species = gas.species_names
        air_species = air.species_names
        MW_O2 = gas.molecular_weights[species.index('O2')]
        MW_N2 = gas.molecular_weights[species.index('N2')]
        MW_CH4 = gas.molecular_weights[species.index('CH4')]
        MW_H2 = gas.molecular_weights[species.index('H2')]
        if self.fuel == 'CH4':
            # CH4 + 2 ( O2 + 3.76 N2) = CO2 + 2 H2O + 3.76(2)N2
            self.LHV = self.LHV_CH4_calc(gas)  # MJ/kg
            self.AF_stoic = 2*(MW_O2+3.76*MW_N2)/MW_CH4

        if self.fuel == 'H2':
            # 2 H2 + ( O2 + 3.76 N2) = 2 H2O + 3.76 N2
            self.LHV = self.LHV_H2_calc(gas)  # MJ/kg
            self.AF_stoic = (MW_O2 + 3.76 * MW_N2) / 2 * MW_H2

        self.mdot_fuel = (self.P_therm/self.LHV)*(1.0/1000.0)  # kg/s
        gas.TPY = 293.15, ct.one_atm, {self.fuel:1.0}  # NTP
        self.vdot_fuel = (self.mdot_fuel/gas.density) * 60000.0  # lnpm
        self.mdot_air = self.AF_stoic * self.mdot_fuel / self.phi
        air.TP = 293.15, ct.one_atm # NTP
        self.vdot_air = (self.mdot_air / air.density) * 60000  # lnpm
        gas_diluent.TPY = 293.15, ct.one_atm, {'N2': 1.0}  # NTP
        self.mdot_N2 = self.vdot_N2 * gas_diluent.density/60000 # kg/s
        self.mdot_total = self.mdot_air + self.mdot_fuel + self.mdot_N2
        air.TP = self.T_air, ct.one_atm
        H_mix = (gas.enthalpy_mass*self.mdot_fuel + air.enthalpy_mass *\
                 self.mdot_air + gas_diluent.enthalpy_mass*self.mdot_N2)/self.mdot_total
        self.enthalpy_inlet = H_mix
        self.Y_comp = {self.fuel: self.mdot_fuel / self.mdot_total,
                       'O2': air.Y[air_species.index('O2')] * self.mdot_air / self.mdot_total,
                       'N2': (air.Y[air_species.index('N2')] * self.mdot_air +self.mdot_N2)/ self.mdot_total,
                       'AR': air.Y[air_species.index('AR')] * self.mdot_air / self.mdot_total}

        gas.HPY = H_mix, ct.one_atm, self.Y_comp
        self.Q_in = H_mix*self.mdot_total
        self.T_in = gas.T

        res1 = ct.Reservoir(gas)
        gas.TPY = 2000,ct.one_atm, self.Y_comp
        react = ct.IdealGasReactor(gas, energy='on')
        vol_comb = (3.14 / 4.0) * (0.215 ** 2.0) * 0.49  # m^3
        react.volume = vol_comb/12.0
        mfc1 = ct.MassFlowController(res1,react,mdot = self.mdot_total)
        res2 = ct.Reservoir(gas)
        mfc2 = ct.MassFlowController(react,res2,mdot = self.mdot_total)
        net = ct.ReactorNet([react])
        dt = 1e-4
        tf = dt
        # net.advance_to_steady_state()
        for i in range(int(5e5)):
            net.advance(tf)
            tf = tf + dt

        self.enthalpy_outlet = react.thermo.enthalpy_mass

        print react.thermo.T
        print "O2=",react.thermo.X[species.index('O2')]
        print "CO2=",react.thermo.X[species.index('CO2')]
        print "CO=", react.thermo.X[species.index('CO')]
        print "Volume settings reactor=",react.volume











