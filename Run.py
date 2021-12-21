import numpy as np
import Settings as st
import BurnerHead as bh
import Combustor as comb
import pickle
import matplotlib.pyplot as plt
import os

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

    def readfiledata(self,case_no):
        case = {1: "0.6", 2: "0.8", 3: "0.8_N2450lnpm"}
        header1 = "_6mmnozzle_TC_60kWphi"
        fileloc2 = "C:/Users/rishikeshsampa/Documents/Research/CombustorGasCompositionandTemperaturemeasurement_November2021/"
        header_loss = "HeatLoss_Port3"
        with open(fileloc2 + header_loss + header1 + case[case_no] + ".dat") as file:
            line = file.readlines()
            for i in range(len(line)):
                if i == 0:
                    continue
                else:
                    ln = line[i].split()
                    T_heat=(float(ln[3]))
        return T_heat

    def piv_data_eval(self):
        """
        Combustor run at constant external dilution but varying equivalence ratio and thermal power
        :return:
        """
        case = {1:"phi06",2:"phi08",3:"phi08_N2"}
        phi = [0.6,0.8,0.8]
        vol_N2=[0.0,0.0,450] #lnpm
        power = 60 #kW
        data_adiabatic={1:{},2:{},3:{}}
        data_nonadiabatic = {1:{},2:{},3:{}}
        for i in range(len(case)):
            T_heater = 273.15+self.readfiledata(i+1)
            settings = st.Settings(self.fuel, power, phi[i], T_heater,vol_N2[i])
            burner = bh.BurnerHead(settings, T_heater)
            combustor = comb.Combustor(settings)
            combustor.adiabatic =True
            data_adiabatic[i+1] = combustor.combustor(i+1)
            # Non-adiabatic cases
            settings2 = st.Settings(self.fuel, power, phi[i], T_heater,vol_N2[i])
            burner2 = bh.BurnerHead(settings2, T_heater)
            combustor2 = comb.Combustor(settings2)
            combustor2.adiabatic =False
            data_nonadiabatic[i + 1] = combustor2.combustor(i + 1)

        os.chdir('C:/Users/rishikeshsampa/Documents/Research/AGNES_proletariat/')
        with open('data_adiabatic.pkl', 'wb') as file:
            pickle.dump(data_adiabatic, file, pickle.HIGHEST_PROTOCOL)

        with open('data_nonadiabatic.pkl', 'wb') as file:
            pickle.dump(data_nonadiabatic, file, pickle.HIGHEST_PROTOCOL)

    def emissplot(self):
        os.chdir('C:/Users/rishikeshsampa/Documents/Research/AGNES_proletariat/')
        with open("data_adiabatic.pkl", 'rb') as f:
            data_adiabatic = pickle.load(f)
        with open("data_nonadiabatic.pkl", 'rb') as f:
            data_nonadiabatic = pickle.load(f)
        cases = data_adiabatic.keys()
        case_name = {1:"$\phi$=0.6",2:"$\phi$=0.8",3:"$\phi$=0.8 N2 dilution"}
        leg = ["NOx", "CO"]
        leg1 = ["CO2", "CH4"]
        leg2 = ["O2", "H2O"]
        X = np.arange(len(data_adiabatic[cases[0]][0]))
        zones = X
        marker_list =["+","*","x"]
        color_list =['k','r','m','g']
        m_count =0
        fig, ax = plt.subplots(dpi=110)
        for i in [0,1,2]:
            data = data_adiabatic[i+1]
            data2 = data_nonadiabatic[i+1]
            ax.scatter(X, data[0],label=case_name[i+1],s = 20,marker=marker_list[m_count],color=color_list[0])
            ax.scatter(X, data2[0], label=case_name[i + 1]+" w/ heat loss",s = 20,marker=marker_list[m_count],color=color_list[i+1])
            m_count += 1
            ax.set_xticks(range(len(zones)))
            ax.set_xticklabels([num+1 for num in zones])
            # ax.bar(X + 0.25, data[1], width=0.25, align='center', tick_label=["%d" % num for num in zones])
            # ax.set_yscale('log')
            # ax.legend(labels=leg)
            ax.set_ylabel("NOx dry mole fraction at 15% O2 (ppm)")
            ax.set_xlabel("Reactor")
            """par1 = ax.twinx()
            par1.bar(X + 0.25, data[1], width=0.25, align='center',
                     tick_label=["%d" % num for num in zones], color="orange", label="CO")
            par1.set_ylabel("CO dry """
            # ask matplotlib for the plotted objects and their labels
            lines, labels = ax.get_legend_handles_labels()
            #lines1, labels1 = par1.get_legend_handles_labels()
            #plt.legend(lines + lines1, labels + labels1, loc=2)
            plt.legend(lines, labels, loc=2, prop={'size': 9}, markerscale=0.85, ncol=3, \
                       bbox_to_anchor=(-0.12, 1.155), columnspacing=0.2)
        ax.set_yscale('log')
        plt.savefig('NOx_comparison.pdf')
        plt.savefig('NOx_comparison.png')

        m_count = 0
        fig2, ax2 = plt.subplots(dpi = 110)
        for i in [0, 1, 2]:
            data = data_adiabatic[i + 1]
            data2 = data_nonadiabatic[i + 1]
            ax2.scatter(X, data[6], label=case_name[i + 1],s = 20,marker=marker_list[m_count],color=color_list[0])
            ax2.scatter(X, data2[6], label=case_name[i + 1]+" w/ heat loss",s = 20,marker=marker_list[m_count],color=color_list[i+1])
            m_count += 1
            ax2.set_xticks(range(len(zones)))
            ax2.set_xticklabels([num+1 for num in zones])
            ax2.set_ylabel("Temperature (K)")
            ax2.set_xlabel("Reactor")
            """par1 = ax.twinx()
            par1.bar(X + 0.25, data[1], width=0.25, align='center',
                     tick_label=["%d" % num for num in zones], color="orange", label="CO")
            par1.set_ylabel("CO dry """
            # ask matplotlib for the plotted objects and their labels
            lines, labels = ax2.get_legend_handles_labels()
            # lines1, labels1 = par1.get_legend_handles_labels()
            # plt.legend(lines + lines1, labels + labels1, loc=2)
            plt.legend(lines, labels, loc=2,prop={'size':9},markerscale=0.85,ncol=3,\
                       bbox_to_anchor=(-0.12,1.155),columnspacing =0.2)
        plt.savefig('Temperature_comparison.pdf')
        plt.savefig('Temperature_comparison.png')

        m_count = 0
        fig3, ax3 = plt.subplots(dpi=110)
        for i in [0, 1, 2]:
            data = data_adiabatic[i + 1]
            data2 = data_nonadiabatic[i + 1]
            ax3.scatter(X, data[1], label=case_name[i + 1], s=20, marker=marker_list[m_count], color=color_list[0])
            ax3.scatter(X, data2[1], label=case_name[i + 1] + " w/ heat loss", s=20, marker=marker_list[m_count],
                        color=color_list[i + 1])
            m_count += 1
            ax3.set_xticks(range(len(zones)))
            ax3.set_xticklabels([num + 1 for num in zones])
            ax3.set_ylabel("CO dry mole fraction at 15% O2 (ppm)")
            ax3.set_xlabel("Reactor")
            ax3.set_yscale('log')
            """par1 = ax.twinx()
            par1.bar(X + 0.25, data[1], width=0.25, align='center',
                     tick_label=["%d" % num for num in zones], color="orange", label="CO")
            par1.set_ylabel("CO dry """
            # ask matplotlib for the plotted objects and their labels
            lines, labels = ax3.get_legend_handles_labels()
            # lines1, labels1 = par1.get_legend_handles_labels()
            # plt.legend(lines + lines1, labels + labels1, loc=2)
            plt.legend(lines, labels, loc=2, prop={'size': 9}, markerscale=0.85, ncol=3, \
                       bbox_to_anchor=(-0.12, 1.155), columnspacing=0.2)
        plt.savefig('CO_comparison.pdf')
        plt.savefig('CO_comparison.png')

        m_count = 0
        fig4, ax4 = plt.subplots(dpi=110)
        for i in [0, 1, 2]:
            data = data_adiabatic[i + 1]
            data2 = data_nonadiabatic[i + 1]
            ax4.scatter(X, data[2], label=case_name[i + 1], s=20, marker=marker_list[m_count], color=color_list[0])
            ax4.scatter(X, data2[2], label=case_name[i + 1] + " w/ heat loss", s=20, marker=marker_list[m_count],
                        color=color_list[i + 1])
            m_count += 1
            ax4.set_xticks(range(len(zones)))
            ax4.set_xticklabels([num + 1 for num in zones])
            ax4.set_ylabel("CO2 dry mole fraction at 15% O2 (ppm)")
            ax4.set_xlabel("Reactor")
            #ax4.set_yscale('log')
            # ask matplotlib for the plotted objects and their labels
            lines, labels = ax4.get_legend_handles_labels()
            # lines1, labels1 = par1.get_legend_handles_labels()
            # plt.legend(lines + lines1, labels + labels1, loc=2)
            plt.legend(lines, labels, loc=2, prop={'size': 9}, markerscale=0.85, ncol=3, \
                       bbox_to_anchor=(-0.12, 1.155), columnspacing=0.2)
        plt.savefig('CO2_comparison.pdf')
        plt.savefig('CO2_comparison.png')

        m_count = 0
        fig5, ax5 = plt.subplots(dpi=110)
        for i in [0, 1, 2]:
            data = data_adiabatic[i + 1]
            data2 = data_nonadiabatic[i + 1]
            ax5.scatter(X, data[3], label=case_name[i + 1], s=20, marker=marker_list[m_count], color=color_list[0])
            ax5.scatter(X, data2[3], label=case_name[i + 1] + " w/ heat loss", s=20, marker=marker_list[m_count],
                        color=color_list[i + 1])
            m_count += 1
            ax5.set_xticks(range(len(zones)))
            ax5.set_xticklabels([num + 1 for num in zones])
            ax5.set_ylabel("CH4 dry mole fraction at 15% O2 (ppm)")
            ax5.set_xlabel("Reactor")
            ax5.set_yscale('log')
            # ask matplotlib for the plotted objects and their labels
            lines, labels = ax5.get_legend_handles_labels()
            # lines1, labels1 = par1.get_legend_handles_labels()
            # plt.legend(lines + lines1, labels + labels1, loc=2)
            plt.legend(lines, labels, loc=2, prop={'size': 9}, markerscale=0.85, ncol=3, \
                       bbox_to_anchor=(-0.12, 1.155), columnspacing=0.2)
        plt.savefig('CH4_comparison.pdf')
        plt.savefig('CH4_comparison.png')

        m_count = 0
        fig6, ax6 = plt.subplots(dpi=115)
        for i in [0, 1, 2]:
            data = data_adiabatic[i + 1]
            data2 = data_nonadiabatic[i + 1]
            ax6.scatter(X, [x*100 for x in data[4]], label=case_name[i + 1], s=20, marker=marker_list[m_count], color=color_list[0])
            ax6.scatter(X, [x*100 for x in data2[4]], label=case_name[i + 1] + " w/ heat loss", s=20, marker=marker_list[m_count],
                        color=color_list[i + 1])
            m_count += 1
            ax6.set_xticks(range(len(zones)))
            ax6.set_xticklabels([num + 1 for num in zones])
            ax6.set_ylabel("O2%")
            ax6.set_xlabel("Reactor")
            #ax6.set_yscale('log')
            # ask matplotlib for the plotted objects and their labels
            lines, labels = ax6.get_legend_handles_labels()
            # lines1, labels1 = par1.get_legend_handles_labels()
            # plt.legend(lines + lines1, labels + labels1, loc=2)
            plt.legend(lines, labels, loc=2, prop={'size': 9}, markerscale=0.85, ncol=3, \
                       bbox_to_anchor=(-0.12, 1.155), columnspacing=0.2)
        plt.savefig('O2_comparison.pdf')
        plt.savefig('O2_comparison.png')
        plt.show()


if __name__ == "__main__":
    run = Run()
    #run.const_dilution()
    #run.piv_data_eval()
    run.emissplot()