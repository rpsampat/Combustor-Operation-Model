import cantera as ct


class Settings:

    def __init__(self, fuel, P_therm, phi, T_heater):
        self.P_therm = P_therm  # kW
        self.fuel = fuel  # 'CH4'/'H2'/'CH4+H2'
        self.LHV = 0
        self.phi = phi
        self.AF_stoic = 0
        self.O2_perc = 0
        self.mdot_fuel = 0  # kg/s
        self.mdot_air = 0 # kg/s
        self.mdot_total = 0 #kg/s
        self.vdot_fuel = 0 # lnpm
        self.vdot_air = 0 # lnpm
        self.T_air = T_heater
        self.T_fuel = 273.15 + 15
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

        self.mdot_fuel = (self.P_therm/self.LHV)*(1/1000.0)  # kg/s
        gas.TPY = 293.15, ct.one_atm, {self.fuel:1.0}  # NTP
        self.vdot_fuel = (self.mdot_fuel/gas.density) * 60000.0  # lnpm
        self.mdot_air = self.AF_stoic * self.mdot_fuel / self.phi
        air.TP = 293.15, ct.one_atm # NTP
        self.vdot_air = (self.mdot_air / air.density) * 60000  # lnpm
        self.mdot_total = self.mdot_air + self.mdot_fuel
        air.TP = self.T_air, ct.one_atm
        H_mix = (gas.enthalpy_mass*self.mdot_fuel + air.enthalpy_mass * self.mdot_air)/self.mdot_total
        gas.HPY = H_mix, ct.one_atm, {self.fuel: self.mdot_fuel/self.mdot_total,
                                      'O2': air.Y[air_species.index('O2')]*self.mdot_air/self.mdot_total,
                                      'N2': air.Y[air_species.index('N2')]*self.mdot_air/self.mdot_total,
                                      'AR': air.Y[air_species.index('AR')]*self.mdot_air/self.mdot_total}









