import cantera as ct
import numpy as np

def premix_freeflame(phi):
    # Simulation parameters
    p = ct.one_atm  # pressure [Pa]
    Tin = 300.0  # unburned gas temperature [K]
    #phi = 1.0
    a = 2.0/phi
    width = 0.1 # m
    loglevel = 1  # amount of diagnostic output (0 to 8)

    # IdealGasMix object used to compute mixture properties, set to the state of the
    # upstream fuel-air mixture
    gas = ct.Solution('gri30.cti')
    #gas = ct.Solution('h2o2.xml')
    reactants = {'CH4':1.0, 'O2':a*1.0, 'N2':a*3.76, 'CO2': 0.0*(a-1.5)}  # premixed gas composition
    #reactants = 'H2:1.1, O2:1, AR:5'  # premixed gas composition
    gas.TPX = Tin, p, reactants

    # Set up flame object
    f = ct.FreeFlame(gas, width=width)
    #f = ct.BurnerFlame(gas, width=width)
    #f=ct.CounterflowPremixedFlame(gas,width=width)
    f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
    #f.show_solution()

    # Solve with mixture-averaged transport model
    f.transport_model = 'Mix'
    f.solve(loglevel=loglevel, auto=True)
    return f

if __name__=='__main__':
    f=premix_freeflame(0.7)
    points=f.flame.n_points
    xloc =  f.grid * 1000.0 # mm
    T=f.T
    g_react = ct.Solution('gri30.cti')
    g=f.gas
    index=g_react.species_index('CO2')
    conc=[]
    dT=[]
    Y_co2 =[]
    dx= 0.03/points
    for i in range(points):
        f.set_gas_state(i)
        molefrac=g.X
        try:
            dT.append((T[i+1]-T[i]))#/dx)
        except:
            pass
        conc.append(molefrac[index])
    import matplotlib.pyplot as plt
    ax=plt.subplot()
    ax.plot(xloc,T, 'r')
    ax1=ax.twinx()
    ax1.plot(xloc,conc)
    plt.show()
