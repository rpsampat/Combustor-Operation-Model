import cantera as ct
import numpy as np
import matplotlib.pyplot as plt


class PathwayAnalysis:
    def __init__(self, combustor_crn):
        # combustor_crn = { zoneid:[PSR1,PSR2..],...}
        zones = combustor_crn.keys()
        gas = combustor_crn[zones[0]][0].thermo # gas object of first reactor of first zone
        self.species = gas.species_names
        self.pathway = {}
        self.combustor_crn = combustor_crn

    def radical_consump(self,PSR):
        # gas in PSR
        g = PSR.thermo
        # number of reactions
        number = g.n_reactions
        k_tot = 0
        reac1 = 'O + CH <=> H + CO'
        reac2 = 'O + CH2 <=> H + HCO'
        reac3 = 'O + CH3 <=> H + CH2O'
        reac4 = 'O + CH3 <=> H + H2 + CO'
        for i in range(number):
            k1 = 0
            k2 = 0
            k = 0
            str = g.reaction_equation(i)
            if str == reac1 or str == reac2 or str == reac3 or str == reac4:
                react = g.reactants(i)
                react = react.split(' ')
                prod = g.products(i)
                prod = prod.split(' ')
                try:
                    # species in product
                    check = prod.index('O')
                    ind = g.kinetics_species_index('O')
                    coeff = g.product_stoich_coeff(ind, i)
                    k = g.net_rates_of_progress[i] * (coeff)
                    k_tot += k
                except:
                    try:
                        # species in reactant
                        check = react.index('O')
                        ind = g.kinetics_species_index('O')
                        coeff = g.reactant_stoich_coeff(ind, i)
                        k = -g.net_rates_of_progress[i] * (coeff)
                        k_tot += k
                    except:
                        continue

        return k_tot

    def prompt_pathway(self,PSR):
        # gas in PSR
        g = PSR.thermo
        # number of reactions
        number = g.n_reactions
        k_tot = 0
        reac1 = 'CH + N2 <=> HCN + N'
        reac2 = 'CH + N2 (+M) <=> HCNN (+M)'
        reac3 = 'CH2 + N2 <=> HCN + NH'
        reac4 = 'H2CN + N <=> N2 + CH2'
        reac5 = 'HCNN + H <=> CH2 + N2'
        reac6 = 'CH2(S) + N2 <=> NH + HCN'
        # reac3='N + OH <=> H + NO'
        for i in range(number):
            k1 = 0
            k2 = 0
            k = 0
            str = g.reaction_equation(i)
            if str == reac1 or str == reac2 or str == reac3 \
                    or str == reac4 or str == reac5 or str == reac6:
                react = g.reactants(i)
                react = react.split(' ')
                prod = g.products(i)
                prod = prod.split(' ')
                try:
                    # species in product
                    check = prod.index('N2')
                    ind = g.kinetics_species_index('N2')
                    coeff = g.product_stoich_coeff(ind, i)
                    k = g.net_rates_of_progress[i] * (coeff)
                    k_tot += k
                except:
                    try:
                        # species in reactant
                        check = react.index('N2')
                        ind = g.kinetics_species_index('N2')
                        coeff = g.reactant_stoich_coeff(ind, i)
                        k = -g.net_rates_of_progress[i] * (coeff)
                        k_tot += k
                    except:
                        continue

        return k_tot

    def reburn_pathway(self,PSR):
        # gas in PSR
        g = PSR.thermo
        # number of reactions
        number = g.n_reactions
        k_tot = 0
        reac1 = 'CH + NO <=> HCN + O'
        reac2 = 'CH + NO <=> H + NCO'
        reac3 = 'CH + NO <=> N + HCO'
        reac4 = 'CH3 + NO <=> HCN + H2O'
        reac5 = 'CH3 + NO <=> H2CN + OH'
        reac6 = 'CH2 + NO <=> OH + HCN'
        reac7 = 'CH2 + NO <=> H + HNCO'
        reac8 = 'CH2 + NO <=> H + HCNO'
        reac9 = 'CH2(S) + NO <=> H + HNCO'
        reac10 = 'CH2(S) + NO <=> OH + HCN'
        reac11 = 'CH2(S) + NO <=> H + HCNO'
        reac12 = 'HCCO + NO <=> HCNO + CO'
        reac13 = 'C + NO <=> CN + O'
        reac14 = 'C + NO <=> CO + N'

        for i in range(number):
            k1 = 0
            k2 = 0
            k = 0
            str = g.reaction_equation(i)
            if str == reac1 or str == reac2 or str == reac3 or str == reac4 or str == reac5 \
                    or str == reac6 or str == reac7 or str == reac8 or str == reac9 \
                    or str == reac10 or str == reac11 or str == reac12 or str == reac13 or str == reac14:
                react = g.reactants(i)
                react = react.split(' ')
                prod = g.products(i)
                prod = prod.split(' ')
                try:
                    # species in product
                    check = prod.index('NO')
                    ind = g.kinetics_species_index('NO')
                    coeff = g.product_stoich_coeff(ind, i)
                    k = g.net_rates_of_progress[i] * (coeff)
                    k_tot += k
                except:
                    try:
                        # species in reactant
                        check = react.index('NO')
                        ind = g.kinetics_species_index('NO')
                        coeff = g.reactant_stoich_coeff(ind, i)
                        k = -g.net_rates_of_progress[i] * (coeff)
                        k_tot += k
                    except:
                        continue

        return k_tot

    def thermal_pathway(self,PSR):
        # gas in PSR
        g = PSR.thermo
        # number of reactions
        number = g.n_reactions
        k_tot = 0
        reac1 = 'N + NO <=> N2 + O'
        # reac2='N + O2 <=> NO + O'
        # reac3='N + OH <=> H + NO'
        for i in range(number):
            k1 = 0
            k2 = 0
            k = 0
            str = g.reaction_equation(i)
            if str == reac1:  # or str == reac2 or str == reac3:
                react = g.reactants(i)
                react = react.split(' ')
                prod = g.products(i)
                prod = prod.split(' ')
                try:
                    # species in product
                    check = prod.index('NO')
                    ind = g.kinetics_species_index('NO')
                    coeff = g.product_stoich_coeff(ind, i)
                    k = g.net_rates_of_progress[i] * (coeff)
                    k_tot += k
                except:
                    try:
                        # species in reactant
                        check = react.index('NO')
                        ind = g.kinetics_species_index('NO')
                        coeff = g.reactant_stoich_coeff(ind, i)
                        k = -g.net_rates_of_progress[i] * (coeff)
                        k_tot += k
                    except:
                        continue

        return k_tot

    def n2o_pathway(self,PSR):
        # gas in PSR
        g = PSR.thermo
        # number of reactions
        number = g.n_reactions
        k_tot = 0
        reac1 = 'N2O (+M) <=> N2 + O (+M)'
        reac2 = 'N2O + O <=> 2 NO'
        for i in range(number):
            k1 = 0
            k2 = 0
            k = 0
            str = g.reaction_equation(i)
            if str == reac1 or str == reac2:
                react = g.reactants(i)
                react = react.split(' ')
                prod = g.products(i)
                prod = prod.split(' ')
                try:
                    # species in product
                    check = prod.index('NO')
                    ind = g.kinetics_species_index('NO')
                    coeff = g.product_stoich_coeff(ind, i)
                    k = g.net_rates_of_progress[i] * (coeff)
                    k_tot += k
                except:
                    try:
                        # species in reactant
                        check = react.index('NO')
                        ind = g.kinetics_species_index('NO')
                        coeff = g.reactant_stoich_coeff(ind, i)
                        k = -g.net_rates_of_progress[i] * (coeff)
                        k_tot += k
                    except:
                        continue

        return k_tot

    def nnh_pathway(self,PSR):
        # gas in PSR
        g = PSR.thermo
        # number of reactions
        number = g.n_reactions
        k_tot = 0
        reac1 = 'NNH <=> H + N2'
        reac2 = 'NNH + O <=> NH + NO'

        for i in range(number):
            k1 = 0
            k2 = 0
            k = 0
            str = g.reaction_equation(i)
            if str == reac1 or str == reac2:
                react = g.reactants(i)
                react = react.split(' ')
                prod = g.products(i)
                prod = prod.split(' ')
                try:
                    # species in product
                    check = prod.index('NO')
                    ind = g.kinetics_species_index('NO')
                    coeff = g.product_stoich_coeff(ind, i)
                    k = g.net_rates_of_progress[i] * (coeff)
                    k_tot += k
                except:
                    try:
                        # species in reactant
                        check = react.index('NO')
                        ind = g.kinetics_species_index('NO')
                        coeff = g.reactant_stoich_coeff(ind, i)
                        k = -g.net_rates_of_progress[i] * (coeff)
                        k_tot += k
                    except:
                        continue

        return k_tot

    def Post(self):
        Queue = []
        diction = {}
        diction['Thermal NOx'] = {}
        diction['Prompt NOx'] = {}
        diction['Reburn NOx'] = {}
        #diction['NOx Destr Rate'] = {}
        #diction['O Alt consump'] = {}
        #diction['NOx Net Rate(bin)'] = {}
        diction['NOx Net Rate'] = {}
        species = self.species
        for zone in self.combustor_crn.keys():
            #diction['O Alt consump'][zone] =[]
            #diction['NOx Destr Rate'][zone] =[]
            diction['NOx Net Rate'][zone] =[]
            #diction['NOx Net Rate(bin)'][zone] =[]
            diction['Thermal NOx'][zone] =[]
            diction['Prompt NOx'][zone] =[]
            diction['Reburn NOx'][zone] =[]
            for reactor in self.combustor_crn[zone]:
                gas = reactor.thermo
                noxprod = gas.creation_rates[species.index('NO')]
                noxdestr = gas.destruction_rates[species.index('NO')] + 1e-20
                try:
                    diction['NOx Destr Rate'][zone].append(noxdestr)
                except:
                    pass
                noxnet = gas.net_production_rates[species.index('NO')] + 1e-20
                try:
                    diction['NOx Net Rate'][zone].append(noxnet)
                except:
                    pass
                # net nox rate
                netrate = 0
                if noxnet > 0:
                    netrate = 1
                elif noxnet < 0:
                    netrate = -1
                try:
                    diction['NOx Net Rate(bin)'][zone].append(netrate)
                except:
                    pass
                # thermal nox
                thermalpath = self.thermal_pathway(reactor)
                tnox = 1e-15
                if thermalpath > 0:
                    tnox = thermalpath
                diction['Thermal NOx'][zone].append(tnox)
                # promptnox
                ppath = self.prompt_pathway(reactor)
                phcn = 1e-20
                if ppath < 0:
                    phcn = -ppath
                diction['Prompt NOx'][zone].append(phcn)
                # nox reburn
                reburnpath = self.reburn_pathway(reactor)
                rbnr = 1e-20
                if reburnpath < 0:  # and noxdestr>0:
                    rbnr = -reburnpath  # /noxdestr
                diction['Reburn NOx'][zone].append(rbnr)
                # O radical consumption
                radicconsump = self.radical_consump(reactor)
                rcr = 1e-20
                if radicconsump < 0:
                    rcr = -radicconsump
                try:
                    diction['O Alt consump'][zone].append(rcr)
                except:
                    pass
        self.pathway = diction

    def plot(self):
        leg= self.pathway.keys()
        spacing = 0.0
        fig, ax = plt.subplots()
        fig2, ax2 = plt.subplots()
        for k in leg:
            dat_samp = []
            X = np.arange(len(self.pathway[k].keys()))
            for zone in self.pathway[k]:
                dat_samp.append(self.pathway[k][zone][-1])
            if k == "NOx Net Rate":
                ax2.bar(X + spacing, dat_samp, width=0.25, align='center',
                       tick_label=["zone %d" % num for num in self.pathway[k].keys()])
            else:
                ax.bar(X + spacing, dat_samp, width = 0.25, align = 'center',
                   tick_label=["zone %d" % num for num in self.pathway[k].keys()])
            spacing = spacing + 0.25

        ax.set_yscale('log')
        leg.remove("NOx Net Rate")
        ax.legend(labels=leg)
        ax.set_ylabel('Rate of NOx production (mol/m^3-s)')
        ax2.set_yscale('log')
        ax2.legend(labels=["NOx Net Rate"])
        ax2.set_ylabel('Rate of NOx production (mol/m^3-s)')
        plt.show()

    def main(self):
        self.Post()
        self.plot()







if __name__ == "__main__":
    loc = os.getcwd()
    print "Reading Data files"
    with open("header.pkl", 'rb') as f:
        header = pickle.load(f)
    with open("data.pkl", 'rb') as f:
        data = pickle.load(f)

    path = loc + '/5000_4827_Static Temperature_Mass fraction of h2o_Mass fraction of ch4_Velocity Dir_doe_reactnum'
    print "Reading Graph"
    with open(path + "/graph2plot.pkl", 'rb') as f:
        graph = pickle.load(f)
    Post(path, graph, header, data)


