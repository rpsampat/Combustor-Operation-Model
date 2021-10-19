import cantera as ct
import math
import numpy as np
from numpy import linalg
from pathway_analysis import PathwayAnalysis


class CRN:

    def __init__(self,settings,Q_loss,wall_area):
        """
        The combustor is divided into zones. Each zone can contain a number of equal reactors connected in series.
        """
        self.zone = {'ingestion': [1], 'flame1': [2], 'flame2': [3], 'PRZ': [4], 'CRZ': [5], 'exhaust': [6]}
        self.crn_def = {1: [1,0.01], 2: [2,0.145], 3: [2,0.145], 4: [3,0.1], 5: [3,0.1], 6: [2,0.5]}
        # {zone id:[number of reactors in zone, zone volume fraction of total combustor volume]}
        self.connect = {1: [2], 2: [3], 3: [4, 5, 6], 4: [2], 5: [2], 6: []}
        self.alpha_mat = {1: [1.0], 2: [1.0], 3: [0.25, 0.25, 0.5], 4: [1.0], 5: [1.0], 6: [1.0]}
        # alpha_mat= { zone id: [fraction of outlet to zone corresponding to list in self.connect.]..}
        """self.crn_def = {1: [1, 0.5], 2: [1, 0.5]}  
        # {zone id:[number of reactors in zone, zone volume fraction of total combustor volume]///}
        self.connect = {1: [2], 2: []}
        self.alpha_mat = {1: [1.0], 2: []}"""
        self.mass_in = settings.mdot_total
        self.T_in = settings.T_in
        self.Y_comp = settings.Y_comp
        self.inlet = {1: self.mass_in}
        self.massflow = self.mass_balance()
        self.outlet = [6]  # list of reactor numbers serving as outlets of the combustor
        self.vol_comb = 0.01778
        self.heat_loss = {4:Q_loss, 6: Q_loss}  # {zoneid:heat loss W/m^2}
        self.wall_area = {4:wall_area/2.0, 6: wall_area/2.0}


    def mass_balance(self):
        n_size = len(self.crn_def.keys())
        coeffmat = np.zeros((n_size,n_size)) # matrix of sized to the number of zones
        for i in range(n_size):
            # col
            for j in range(n_size):
                if j==i:
                    coeffmat[j][i] = 1.0
                elif (j+1) in self.connect[i+1]:
                    ind = self.connect[i+1].index(j+1)
                    coeffmat[j][i] = -1.0*self.alpha_mat[i+1][ind]

        rhsvect = np.zeros(n_size)
        for k in self.inlet.keys():
            rhsvect[k-1] = self.inlet[k]
        x = linalg.solve(coeffmat, rhsvect)
        return x

    def emiss(self,psr):
        NOx_out = []
        N2O_out = []
        CO_out = []
        NH_out = []
        CO2_out = []
        CH2O = []
        T_out = []
        # Correcting CRN emissions
        Nox_sum = 0
        N2O_sum = 0
        NH_sum = 0
        co_sum = 0
        co2_sum = 0
        T_sum = 0
        m_sum = 0
        N_corr_sum = 0
        N_sum = 0
        gas = psr[0].thermo
        species = gas.species_names
        R=8.314
        perc_corr=0.15
        for n in psr:
            P = n.thermo.P
            T = n.thermo.T
            V = n.volume
            D = n.thermo.density
            N_wet = P * V / (R * T)
            H2O = n.thermo.X[species.index('H2O')]
            NO = n.thermo.X[species.index('NO')]
            N2O = n.thermo.X[species.index('N2O')]
            NO2 = n.thermo.X[species.index('NO2')]
            O2 = n.thermo.X[species.index('O2')]
            OH = n.thermo.X[species.index('OH')]
            CO = n.thermo.X[species.index('CO')]
            NH = n.thermo.X[species.index('NH')]
            CO2 = n.thermo.X[species.index('CO2')]
            O2 = n.thermo.X[species.index('O2')]
            CH4 = n.thermo.X[species.index('CH4')]
            N_dry = N_wet * (1 - H2O)
            NO_dry = NO * N_wet / N_dry
            N2O_dry = N2O * N_wet / N_dry
            NO2_dry = NO2 * N_wet / N_dry
            OH_dry = OH * N_wet / N_dry
            CO_dry = CO * N_wet / N_dry
            NH_dry = NH * N_wet / N_dry
            CO2_dry = CO2 * N_wet / N_dry
            O2_dry = O2 * N_wet / N_dry
            CH4_dry = CH4 * N_wet / N_dry
            # O2 correction
            N_O2 = O2 * N_wet
            #print "O2=", N_O2
            O2_corr = (perc_corr / (1 - perc_corr)) * (N_dry - N_O2)
            #print "O2corr=", O2_corr
            N_corr = N_dry - N_O2 + O2_corr
            #print "O2=", N_dry
            N_corr_sum += N_corr
            N_sum += N_wet
            NO_corr = NO_dry * N_dry / N_corr
            N2O_corr = N2O_dry * N_dry / N_corr
            NO2_corr = NO2_dry * N_dry / N_corr
            OH_corr = OH_dry * N_dry / N_corr
            CO_corr = CO_dry * N_dry / N_corr
            NH_corr = NH_dry * N_dry / N_corr
            CO2_corr = CO2_dry * N_dry / N_corr
            CH4_corr = CH4_dry * N_dry / N_corr
            Nox_sum += (NO_corr + N2O_corr + NO2_corr) * N_corr
            N2O_sum += (N2O_corr) * N_corr
            co_sum += (CO_corr) * N_corr
            NH_sum += (NH_corr) * N_corr
            co2_sum += CO2_corr * N_corr
            T_sum += T * D * V
            m_sum += D * V
            # dat_plt['OH'][n]=OH_corr*1e6
        NOx_out.append(Nox_sum * 1e6 / N_corr_sum)
        CO_out.append(co_sum * 1e6 / N_corr_sum)
        CO2_out.append(co2_sum / N_sum)
        T_out.append(T_sum / m_sum)
        return NOx_out,CO_out

    def combustor(self):
        gas = ct.Solution('gri30.cti')
        air = ct.Solution('air.cti')
        ener = 'on'
        vol_total = self.vol_comb
        v_fact=1000.0
        valve=[]
        mfc_rec =[]
        #Reservoirs
        gas.TPX = 300, ct.one_atm, self.Y_comp  # environment gas state definition
        res_env = ct.Reservoir(gas)  # environment for heat loss at wall
        gas.TPX = self.T_in, ct.one_atm, self.Y_comp  # inlet gas state definition
        res_inlet = ct.Reservoir(gas)
        res_outlet = ct.Reservoir(gas)

        psr_test = ct.IdealGasReactor(gas, energy=ener)
        net_test = ct.ReactorNet([psr_test])
        net_test.advance_to_steady_state()
        gas = psr_test.thermo


        # Generating dictionary of PSR object
        combustor_crn = {} # {zone:[PSR1,PSR2...],...}
        reactor_net = []
        for k in self.crn_def.keys():
            n_psr = len(self.crn_def[k])
            for j in range(self.crn_def[k][0]):
                psr = ct.IdealGasReactor(gas, energy=ener)
                psr.volume = vol_total*self.crn_def[k][1]/self.crn_def[k][0]
                if k in self.heat_loss.keys():
                    wall_area = self.wall_area[k]/n_psr
                    wall = ct.Wall(psr,res_env,A=wall_area,U=self.heat_loss[k])
                try:
                    combustor_crn[k].append(psr)
                except:
                    combustor_crn[k]=[psr]
                reactor_net.append(psr)

        # Generating mass flow controllers
        # Note: This implementation does not allow an inlet zone to accept an internal recirculation
        for i in self.connect.keys():
            # connect last reactor of current zone 'i' to first reactor of neighbouring zone
            psr0 = combustor_crn[i][-1]
            for j in range(len(self.connect[i])):
                mflow = self.massflow[i - 1] * self.alpha_mat[i][j]
                k = self.connect[i][j]
                n_psr = len(combustor_crn[k])
                psr1 = combustor_crn[k][0]
                for m in range(n_psr):
                    # connecting psr within a zone with mfcs of equal mass flow in series
                    psr2 = combustor_crn[k][m]
                    if m == 0:
                        mfc = ct.MassFlowController(psr0,psr1,mdot=mflow)
                        mfc_rec.append(mfc)
                        """val = ct.Valve(psr0,psr1)
                        kv = mflow / v_fact
                        val.set_valve_coeff(kv)"""
                    else:
                        mfc = ct.MassFlowController(psr1, psr2, mdot=mflow)
                        mfc_rec.append(mfc)
                        """val = ct.Valve(psr1, psr2)
                        kv = mflow / v_fact
                        val.set_valve_coeff(kv)"""
                    psr1 = psr2
        # inlets
        for inlet in self.inlet.keys():
            psr0 = res_inlet
            n_psr = len(combustor_crn[inlet])  # number of psrs in zone
            psr1 = combustor_crn[inlet][0]   # first reactor of inlet zone
            mflow = self.inlet[inlet]
            for m in range(n_psr):
                # connecting psr within a zone with mfcs of equal mass flow in series
                psr2 = combustor_crn[inlet][m]
                if m == 0:
                    mfc = ct.MassFlowController(psr0, psr1, mdot=mflow)
                    mfc_rec.append(mfc)
                    """val = ct.Valve(psr0, psr1)
                    kv = mflow / v_fact
                    val.set_valve_coeff(kv)"""
                else:
                    mfc = ct.MassFlowController(psr1, psr2, mdot=mflow)
                    mfc_rec.append(mfc)
                    """val = ct.Valve(psr1, psr2)
                    kv = mflow / v_fact
                    val.set_valve_coeff(kv)"""
                psr1 = psr2
        # outlets
        for outlet in self.outlet:
            # connecting last reactor of zone to outlet reservoir
            mfc = ct.MassFlowController(combustor_crn[outlet][-1], res_outlet, mdot=self.massflow[outlet-1])
            mfc_rec.append(mfc)
            """val = ct.Valve(combustor_crn[outlet][-1], res_outlet)
            kv = self.massflow[outlet-1] / v_fact
            val.set_valve_coeff(kv)"""

        net = ct.ReactorNet(reactor_net)
        print combustor_crn
        print "MFCs=",len(mfc_rec)
        dt = 1e-4
        tf = dt
        #net.advance_to_steady_state()
        for i in range(50000):
            net.advance(tf)
            tf = tf + dt

        nox,co = self.emiss([reactor_net[5]])
        print reactor_net[5].thermo.T
        return reactor_net[5].thermo.T,nox[0],co[0], combustor_crn

    def manualcrn(self):
        """
        Manually generated CRN with cantera objects. Typically 2 PSRs to verify if code works.
        :return:
        """
        # defining network
        gas = ct.Solution('gri30.cti')
        #gas.TPY = 294, ct.one_atm, {'CH4': 0.16, 'O2': (1 - 0.16) * 0.23, 'N2': (1 - 0.16) * 0.77}
        # gas = ct.Solution('h2o2.cti')
        gas.TPX = 1000.0, ct.one_atm, {'CH4': 0.16, 'O2': (1 - 0.16) * 0.23, 'N2': (1 - 0.16) * 0.77}
        res_inlet = ct.Reservoir(gas)
        res_outlet = ct.Reservoir(gas)
        psr_test = ct.IdealGasReactor(gas, energy='on')
        net_test = ct.ReactorNet([psr_test])
        net_test.advance_to_steady_state()
        gas = psr_test.thermo
        psr1 = ct.IdealGasReactor(gas, energy='on')
        psr1.volume = self.vol_comb/2.0
        #gas.TPY = 1500, ct.one_atm, {'CH4': 0.16, 'O2': (1 - 0.16) * 0.23, 'N2': (1 - 0.16) * 0.77}
        psr2 = ct.IdealGasReactor(gas, energy='on')
        psr2.volume = self.vol_comb/2.0
        #gas.TPY = 600, ct.one_atm, {'CH4': 0.16, 'O2': (1 - 0.16) * 0.23, 'N2': (1 - 0.16) * 0.77}
        mfc1 = ct.MassFlowController(psr1, psr2, mdot=self.mass_in)
        mfc2 = ct.MassFlowController(res_inlet, psr1, mdot=self.mass_in)
        mfc3 = ct.MassFlowController(psr2, res_outlet, mdot=self.mass_in)
        rn = ct.ReactorNet([psr1,psr2])
        dt = 1e-4
        tf = dt
        # net.advance_to_steady_state()
        for i in range(50000):
            rn.advance(tf)
            tf = tf + dt

        nox, co = self.emiss([psr2])
        print psr2.thermo.T
        return psr2.thermo.T, nox[0], co[0]

    def main(self):
        t,n,c, crn = self.combustor()
        PA = PathwayAnalysis(crn)
        PA.main()

        return t,n,c


if __name__=="__main__":
    import Settings as settings
    T=[]
    nox=[]
    co=[]
    set = settings.Settings()
    set.main()
    crn = CRN(set)
    t, n,c = crn.manualcrn()
    T.append(t)
    nox.append(n)
    co.append(c)
    t,n,c = crn.combustor()
    T.append(t)
    nox.append(n)
    co.append(c)
    print T
    print nox
    print co