import matplotlib.pyplot as plt
import os


class HeatLossModel:
    def __init__(self):
        self.T_port3=0
        self.T_port5 =0
        self.T_flange =0
        self.T_exh_flange = 0
        self.case={1:"0.6",2:"0.8",3:"0.8_N2450lnpm"}
        self.wall_xloc = [0.0,140,240,490]
        self.T_inter =0
        self.epsi_air = 0.009  # air emissivity
        self.epsi_wall = 0.79  # combustor wall emissivity
        self.sigma = 5.6704e-8  # stefan-boltzman constant
        self.convective_coeff = 45
        self.T_film = 50+273.15 # K
        self.Qloss = 0


    def readwall_Temp(self,case_no):
        fileloc1 ="C:/Users/rishikeshsampa/Documents/Research/CombustorGasCompositionandTemperaturemeasurement_November2021/Temperature measurments/Corrected/"
        header1 = "_6mmnozzle_TC_60kWphi"
        for j in [3,5]:
            with open(fileloc1+"Port"+str(j)+header1+self.case[case_no]+".dat") as file:
                line = file.readlines()
                for i in range(len(line)):
                    if i==1:
                        ln = line[i].split()
                        if j==3:
                            self.T_port3 = float(ln[3])
                        else:
                            self.T_port5 = float(ln[3])
                    if i==4:
                        ln = line[i].split()
                        T_gas3 = float(ln[3])
        fileloc2 = "C:/Users/rishikeshsampa/Documents/Research/CombustorGasCompositionandTemperaturemeasurement_November2021/"
        header_loss = "HeatLoss_Port3"
        with open(fileloc2+header_loss+header1+self.case[case_no]+".dat") as file:
            line = file.readlines()
            for i in range(len(line)):
                if i==0:
                    continue
                else:
                    ln = line[i].split()
                    self.T_flange = float(ln[0])
                    T_exh = float(ln[1])
        self.T_exh_flange = ((self.T_port3+273.15)* (T_exh+273.15)/(T_gas3+273.15))-273.15
        return 0

    def temp_interpolate(self,xloc):
        T_base = [self.T_flange,self.T_port3,self.T_port5,self.T_exh_flange]
        for i in range(len(self.wall_xloc)):
            if xloc<self.wall_xloc[i]:
                break
        self.T_inter = T_base[i-1]+(T_base[i]-T_base[i-1])*(xloc-self.wall_xloc[i-1])/(self.wall_xloc[i]-self.wall_xloc[i-1])
        return self.T_inter

    def heatloss(self):
        T_surr = 20+273.15
        T_wall = self.T_inter + 273.15
        Qrad_out = self.sigma*self.epsi_wall*(T_wall**4.0)
        Qrad_in = self.epsi_wall*self.sigma*self.epsi_air*(T_surr**4.0)
        Qconv = self.convective_coeff*(T_wall-self.T_film)
        self.Qloss = Qrad_out-Qrad_in+Qconv

    def temperature_plots(self):
        fileloc1 = "C:/Users/rishikeshsampa/Documents/Research/CombustorGasCompositionandTemperaturemeasurement_November2021/Temperature measurments/Corrected/"
        header1 = "_6mmnozzle_TC_60kWphi"
        case_name = {1: "$\phi$=0.6", 2: "$\phi$=0.8", 3: "$\phi$=0.8 N2 dilution"}
        marker_list = ["+", "*", "x"]
        for j in [3, 5]:
            fig, ax = plt.subplots(dpi=110)

            for case_num in range(3):
                pos = []
                temp_corr = []
                with open(fileloc1 + "Port" + str(j) + header1 + self.case[case_num+1] + ".dat") as file:
                    line = file.readlines()
                    for i in range(len(line)):
                        if i > 0:
                            ln = line[i].split()
                            pos.append(103.5-float(ln[0]))
                            temp_corr.append(float(ln[3])+273.15)
                ax.scatter(pos, temp_corr, label=case_name[case_num + 1], s=25,marker=marker_list[case_num])
                lines, labels = ax.get_legend_handles_labels()
                # lines1, labels1 = par1.get_legend_handles_labels()
                # plt.legend(lines + lines1, labels + labels1, loc=2)
                plt.legend(lines, labels, loc='best', prop={'size': 9}, markerscale=0.95, ncol=1)

            ax.set_ylabel("Temperature (K)")
            ax.set_xlabel("Radial distance from centerline (mm)")
            plt.savefig('Temperature_Port'+str(j)+'.pdf')
            plt.savefig('Temperature_Port'+str(j)+'.png')
            plt.show()

    def main(self, xloc,case_num):
        self.readwall_Temp(case_num)
        T_res = self.temp_interpolate(xloc)
        self.heatloss()



if __name__=="__main__":
    comb_hl = HeatLossModel()
    #comb_hl.main(55,2)
    #print 1
    comb_hl.temperature_plots()



