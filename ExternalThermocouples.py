from nptdms import TdmsFile
import os
import math
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits import axisartist
import Thermocouple as TC

class ExtThermocouple:
    def __init__(self):
        # 'TC9'=='Flang Temp', 'TC11'=='Frame Temp'
        self.quant = ["Flang Temp","Frame Temp","TC15","TC14","TC6"]#,"Avg Exhaust Temp"]
        #self.quant = ["","","""TC15", "TC14", "TC6"]
        self.filename = ["2021_08_13_16_47_32Thermocouplesfromwall_0mm_40kw_gasanalyser.tdms",
                         "2021_08_13_16_50_08Thermocouplesfromwall_1mm_40kw_gasanalyser.tdms",
                         "2021_08_13_16_52_46Thermocouplesfromwall_2mm_40kw_gasanalyser.tdms",
                         "2021_08_13_16_56_16Thermocouplesfromwall_3mm_40kw_gasanalyser.tdms",
                         "2021_08_13_17_00_54Thermocouplesfromwall_4mm_40kw_gasanalyser.tdms",
                         "2021_08_13_17_05_56Thermocouplesfromwall_5mm_40kw_gasanalyser.tdms"]
        self.air = ct.Solution('air.cti')
    def read_tdms(self,file):
        tdms_file = TdmsFile.read(file)
        group = tdms_file._channel_data
        return group

    def data_extract(self, header, group):
        str_add = "Raw"
        for i in group.keys():
            res = header in i
            if res == True and str_add in i:
                break

        return i

    def main(self):
        group = self.read_tdms(self.filename[0])
        key =[]
        for q in self.quant:
            key.append(self.data_extract(q,group))

        fig, ax = plt.subplots()
        for j in key:
            ax.plot(group[j].data)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Temperature (C)")
        ax.legend(self.quant)
        plt.grid(True)
        plt.show()

    def main_gradient(self):
        val ={}
        val_single = {}
        r_loc = [0, 1, 2, 3, 4, 5] # in mm
        h_fac = [0.50, 0.50, 1.0, 1.0, 1.0]
        r_fac = [0.05, 0.05, 0.95, 0.95, 0.95, 0.95]
        h=[]
        tc = TC.Thermocouple()
        tc.dia_tc = 0.0025  # thermocouple diameter
        tc.mdot_cool = 3000  # air flow rate (lnpm)
        tc.T_air = 100 + 273.15  # K
        tc.area_conv = tc.conv_area_calc_ext()
        area_cond = (math.pi / 4.0) * tc.dia_tc ** 2.0
        T_wall = {}
        for q in range(len(self.quant)):
            ind = 0
            val[q] ={}
            for f in self.filename:
                group = self.read_tdms(f)
                j = self.data_extract(self.quant[q], group)
                dat = np.array(group[j].data)
                val[q][r_loc[ind]]= dat
                len_dat = len(dat)
                try:
                    val_single[q].append(np.mean(dat[int(len_dat/2):len_dat]))
                except:
                    val_single[q]=[np.mean(dat[int(len_dat /2):len_dat])]

                if ind == 0:
                    roots = tc.heat_balance_Twall(val_single[q][0] + 273.15, area_cond=area_cond, h_fac=h_fac[q]*r_fac[ind])
                    #print roots
                    #print T_wall
                    #print q
                    T_wall[q] = roots[0]  # K
                else:
                    roots = tc.heat_balance_fluid(T_wall[q] , val_single[q][ind] + 273.15, area_cond=area_cond, h_fac=h_fac[q]*r_fac[ind])
                #print T_wall[q]
                #print val_single[q][ind]
                try:
                    val_single[q][ind] = roots[0] # K
                    #val_single[q][ind] = val_single[q][ind] + 273.15  # K
                except:
                    val_single[q][ind] =  val_single[q][ind] + 273.15 # K
                ind += 1
            dtdy_0 = (val_single[q][1]-val_single[q][0])/((r_loc[1]-r_loc[0])/1000.0)
            T_ref = (val_single[q][0]+val_single[q][1]+val_single[q][2])/3.0  # K
            self.air.TP = T_ref, ct.one_atm
            k = self.air.thermal_conductivity
            h.append(np.abs(k*dtdy_0/(T_ref-T_wall[q])))
        fig, ax = plt.subplots()
        for j in val[0]:
            ax.plot(val[0][j], label=j)  # TODO: x axis time index from data
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Temperature (C)")
        ax.legend(["0 mm", "1 mm", "2 mm", "3 mm", "4 mm", "5 mm"])
        plt.grid(True)

        fig2, ax2 = plt.subplots()
        for j in val_single:
            ax2.plot(r_loc,np.array(val_single[j]) - 273.15)  # TODO: x axis time index from data
        ax2.set_xlabel("Radial distance (mm)")
        ax2.set_ylabel("Temperature (C)")
        ax2.legend(["TC1", "TC2", "TC3", "TC4", "TC5"])
        plt.grid(True)

        fig3, ax3 = plt.subplots()
        ax3.plot(range(len(self.quant)), h)  # TODO: x axis time index from data
        ax3.set_xlabel("Axial distance (mm)")
        ax3.set_ylabel("Convective heat transfer coefficient (W/K-m^2)")
        plt.grid(True)

        fig4, ax4 = plt.subplots()
        for j in val[3]:
            ax4.plot(val[3][j], label=j)  # TODO: x axis time index from data
        ax4.set_xlabel("Time (s)")
        ax4.set_ylabel("Temperature (C)")
        ax4.legend(["0 mm", "1 mm", "2 mm", "3 mm", "4 mm", "5 mm"])
        plt.grid(True)
        plt.show()


if __name__=="__main__":
    exttherm = ExtThermocouple()
    #exttherm.main()
    exttherm.main_gradient()