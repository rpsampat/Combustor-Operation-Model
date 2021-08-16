import cantera as ct
import math

class CRN:

    def __init__(self):

    def emiss(psr):
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
        species = ['H2', 'H', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2', 'C', 'CH', \
                   'CH2', 'CH2(S)', 'CH3', 'CH4', 'CO', 'CO2', 'HCO', 'CH2O', \
                   'CH2OH', 'CH3O', 'CH3OH', 'C2H', 'C2H2', 'C2H3', 'C2H4', 'C2H5', \
                   'C2H6', 'HCCO', 'CH2CO', 'HCCOH', 'N', 'NH', 'NH2', 'NH3', 'NNH', \
                   'NO', 'NO2', 'N2O', 'HNO', 'CN', 'HCN', 'H2CN', 'HCNN', 'HCNO', 'HOCN', \
                   'HNCO', 'NCO', 'N2', 'AR', 'C3H7', 'C3H8', 'CH2CHO', 'CH3CHO']
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

    def connectivity(self):
        crn = {1:[2], 2:[3], 3:[4,5,6],4:[2], 5:[2]}
        inlet = {'fuel': [1], 'air':[1]}
        outlet = {6}
        heat_loss={4,5,6}
        zone = {'ingestion':[1], 'flame1':[2], 'flame2':[3], 'CRZ':[5], 'PRZ':[4], 'exhaust':[6]}
        return crn, inlet, outlet, heat_loss
    def combustor(m_fuel,m_air, recirc_ratio, hl1, hl2):
        gas = ct.Solution('gri30.cti')
        air = ct.Solution('air.cti')
        ener = 'on'
        vol = 0.00039519
        m1 = m_fuel#2.56e-5
        m2 = m_air#0.00057425
        m3 = m1 + m2
        fact=100
        valve=[]
        #Fuel res
        gas.TPY = 300, ct.one_atm, {'CH4': 1}
        res_f=ct.Reservoir(gas)
        #Air res
        gas.TPX = 673.15, ct.one_atm, air.X
        res_a = ct.Reservoir(gas)
        #exhaust res
        gas.TPX = 300, ct.one_atm, air.X
        res_ex = ct.Reservoir(gas)
        #PSR
        gas.TPX = 2500, ct.one_atm, air.X
        psr1=ct.IdealGasReactor(gas, energy=ener)
        v1 = 5.56042E-05*0.9
        #r=gas.destruction_rates
        #print r
        psr1.volume=v1
        psr9 = ct.IdealGasReactor(gas, energy=ener)
        v9 = 5.56042E-05*0.1
        psr1.volume = v9
        psr2 = ct.IdealGasReactor(gas, energy=ener)
        v2 = 5.56042E-05
        psr2.volume=v2
        psr3 = ct.IdealGasReactor(gas, energy=ener)
        v3 = 5.56042E-05
        psr3.volume=v3
        psr4 = ct.IdealGasReactor(gas, energy=ener)
        v4 = 2.05217E-05
        psr4.volume=v4
        w4 = ct.Wall(psr4, res_ex, Q=hl2/3)
        psr5 = ct.IdealGasReactor(gas, energy=ener)
        v5 = 2.05217E-05
        psr5.volume=v5
        w5 = ct.Wall(psr5, res_ex, Q=hl2/3)
        psr6 = ct.IdealGasReactor(gas, energy=ener)
        v6 = 2.05217E-05
        psr6.volume=v6
        w6 = ct.Wall(psr6, res_ex, Q=hl2/3)
        psr7 = ct.IdealGasReactor(gas, energy=ener)
        v7 = 0.000166813/2
        psr7.volume=v7
        w7 = ct.Wall(psr7, res_ex, Q=hl1/2)
        psr8 = ct.IdealGasReactor(gas, energy=ener)
        v8 = 0.000166813/2
        psr8.volume = v8
        w8 = ct.Wall(psr8, res_ex, Q=hl1/2)
        #Recirculation flow
        m_recirc = m3 * recirc_ratio
        m47 = m_recirc * 0.6
        m37 = m_recirc * 0.3
        m28 = m_recirc * 0.1
        m78 = m47 + m37
        m8out = m78 + m28
        m81 = 0.99*m8out
        m82 = 0.01*m8out
        m9in = m1
        m91 = 0.99* m9in
        m92 = 0.01* m9in
        m1in = m2+m91+m81
        m12 = m1in
        m2in = m12+m82
        m23 = m2in-m28
        m3in = m23
        m34 = m3in-m37
        m4in = m34
        m45 = m4in- m47
        m56=m45
        mfc01=ct.MassFlowController(res_f,psr9,mdot=m1)
        val01 = ct.Valve(res_f, psr9)
        kv= m1/fact
        val01.set_valve_coeff(kv)
        mfc91 = ct.MassFlowController(psr9, psr1, mdot=m91)
        val91 = ct.Valve(psr9, psr1)
        kv = m91 / fact
        val91.set_valve_coeff(kv)
        mfc92 = ct.MassFlowController(psr9, psr2, mdot=m92)
        val92 = ct.Valve(psr9, psr2)
        kv = m92 / fact
        val92.set_valve_coeff(kv)
        mfc11 = ct.MassFlowController(res_a, psr1, mdot=m2)
        val11 = ct.Valve(res_a, psr1)
        kv = m2 / fact
        val11.set_valve_coeff(kv)
        mfc12 = ct.MassFlowController(psr1, psr2, mdot=m12)
        val12 = ct.Valve(psr1,psr2)
        kv = m12 / fact
        val12.set_valve_coeff(kv)
        mfc23 = ct.MassFlowController(psr2, psr3, mdot=m23)
        val23 = ct.Valve(psr2, psr3)
        kv = m23 / fact
        val23.set_valve_coeff(kv)
        mfc34 = ct.MassFlowController(psr3, psr4, mdot=m34)
        val34 = ct.Valve(psr3, psr4)
        kv = m34 / fact
        val34.set_valve_coeff(kv)
        mfc45 = ct.MassFlowController(psr4, psr5, mdot=m45)
        val45 = ct.Valve(psr4, psr5)
        kv = m45 / fact
        val45.set_valve_coeff(kv)
        mfc56 = ct.MassFlowController(psr5, psr6, mdot=m56)
        val56 = ct.Valve(psr5, psr6)
        kv = m56 / fact
        val56.set_valve_coeff(kv)
        mfc6ex = ct.MassFlowController(psr6, res_ex, mdot=m3)
        val6ex = ct.Valve(psr6, res_ex)
        kv = m3 / fact
        val6ex.set_valve_coeff(kv)
        #Recirculating flows
        mfc47 = ct.MassFlowController(psr4, psr7, mdot=m47)
        val47 = ct.Valve(psr4, psr7)
        kv = m47 / fact
        val47.set_valve_coeff(kv)
        mfc37 = ct.MassFlowController(psr3, psr7, mdot=m37)
        val37 = ct.Valve(psr3, psr7)
        kv = m37 / fact
        val37.set_valve_coeff(kv)
        mfc28 = ct.MassFlowController(psr2, psr8, mdot=m28)
        val28 = ct.Valve(psr2, psr8)
        kv = m28 / fact
        val28.set_valve_coeff(kv)
        mfc78 = ct.MassFlowController(psr7, psr8, mdot=m78)
        val78 = ct.Valve(psr7, psr8)
        kv = m78 / fact
        val78.set_valve_coeff(kv)
        mfc81 = ct.MassFlowController(psr8, psr1, mdot=m81)
        val81 = ct.Valve(psr8, psr9)
        kv = m81 / fact
        val81.set_valve_coeff(kv)
        mfc82 = ct.MassFlowController(psr8, psr2, mdot=m82)
        val82 = ct.Valve(psr8, psr2)
        kv = m82 / fact
        val82.set_valve_coeff(kv)

        PSR=[psr1,psr2,psr3,psr4,psr5,psr6,psr7,psr8,psr9]
        gas.TPX = 300.0, ct.one_atm, 'H:1.0'
        net = ct.ReactorNet(PSR)
        dt=0.0001
        tf=dt

        valve=[val01,val11,val12,val23,val34,val45,val56,val6ex,val47,val37,val28,val78,val81]
        for i in range(50000):
            e_max = 0
            xin = []
            xfin = []
            for p in PSR:
                xin += list(p.thermo.Y)
            net.advance(tf)
            for p in PSR:
                xfin += list(p.thermo.Y)
            tf = tf + dt
            func = []
            for n in range(len(xin)):
                func.append(xfin[n] - xin[n])
            e = []
            e_div = max(xin)
            for ind in range(len(func)):
                # error_cum += (abs(func[ind]) / (xin[ind] + tol)) ** 2
                e.append((abs(func[ind]) / (e_div)))
            e_rel = max(e)
            #print "Error total=", sum(e)
            index_val = e.index(e_rel)
            #print "Index=", index_val
            #print "Y val=", xfin[index_val]
            if e_rel > e_max:
                e_max = e_rel
            #print "Error=", e_max
            #if e_max <1e-06:
                #for v in valve:
                    #v.set_valve_coeff(0)
            if e_max < 3e-09:
                break

        nox,co=emiss([PSR[5]])
        print PSR[5].thermo.T
        return PSR[5].thermo.T,nox[0],co[0]

if __name__=="__main__":
    mf = 2.56e-5
    ma = [0.000493763,0.000574205,0.000642292,0.000725524,0.000829238,0.000909854,0.000940206]
    recirc=[2.269152216,2.294833333,2.195823353,1.944513981,1.918549708,1.89971123,1.914171843]
    hl1=[456.6363693,418.4599281,387.7261415,362.757005,320.1807849,297.3384848,284.6751501]
    hl2=[249.4737448,228.6168433,211.8260808,198.1847145,174.9240857,162.444672,155.5263236]
    T=[]
    nox=[]
    co=[]
    for i in range(7):
        t,n,c=combustor(mf,ma[i],recirc[i],hl1[i],hl2[i])
        T.append(t)
        nox.append(n)
        co.append(c)
    print T
    print nox
    print co