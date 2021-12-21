"""
Read data files that contain CRN data generated based on PIV data. THese files were originally generated
using a Matlab script
"""
import numpy as np
import math
import os
import Combustor_CRN as CombCRN

class CRNPIV_data:
    def __init__(self):
        self.facelist={}
        self.nodes={}
        self.reactors={}
        self.velflux={}
        self.commonfaces={}
        self.zeroadjustface={}
        self.connect={}
        self.alpha_mat = {}
        self.wall_face ={}
        self.reactor_wall={}
        self.outlet = {14:[]}
        self.volume =0
        self.Dia = 0.215
        self.vol_comb_sector = (3.14 / 4.0) * (self.Dia ** 2.0) * 0.49 *(1/12.0) # m^3

    def readata(self):
        # reading nodes
        with open("crn_nodes.dat") as file:
            line = file.readlines()
            for i in range(len(line)):
                if i==0:
                    continue
                else:
                    ln = line[i].split()
                    self.nodes[float(ln[0])]=ln[1:]
        # reading reactor volumes
        ang_sweep = (2.0*math.pi/12.0) # m
        ang_sweep_nozzle = 0.00667#math.pi*0.00667/2.0#m
        near_jet = [1,2,3,4]
        num_react = 1 # number of reactors in zone
        outlet_react_area = (0.49 - 0.221) * (0.103)*ang_sweep
        vol_jet = 0
        with open("crnvolume.dat") as file:
            line = file.readlines()
            for i in range(len(line)):
                if i==0:
                    continue
                else:
                    ln = line[i].split()
                    reactor_id = int(float(ln[0]))
                    sweep_rad = float(ln[2])*1e-3
                    if reactor_id in near_jet:
                        vol_curr = float(ln[1])*1e-6*ang_sweep_nozzle
                        vol_jet = vol_jet + vol_curr
                    else:
                        num_react = 1
                        vol_curr = float(ln[1]) * 1e-6 * ang_sweep * sweep_rad
                    self.volume = self.volume + vol_curr
                    self.reactors[reactor_id] = [num_react, vol_curr, float(ln[3])] # [number of reactors in zone, volume of zone, xdistance from burner head]
                    self.connect[reactor_id] = []
        vol_eff = self.vol_comb_sector-vol_jet
        for i in self.reactors:
            if i in near_jet:
                self.reactors[i][1] = self.reactors[i][1] / self.vol_comb_sector
            else:
                self.reactors[i][1] = (vol_eff*self.reactors[i][1]/self.volume)/self.vol_comb_sector
        self.volume = self.volume+vol_jet
        self.outlet = {max(self.reactors.keys()):[]}
        # reading velocity flux across faces
        with open("interface_velocityflux.dat") as file:
            line = file.readlines()
            for i in range(len(line)):
                if i == 0:
                    continue
                else:
                    ln = line[i].split()
                    self.velflux[ln[0]+ln[1]] = float(ln[2])
        # reading wall faces
        with open("walledges.dat") as file:
            line = file.readlines()
            for i in range(len(line)):
                if i==0:
                    continue
                else:
                    ln = line[i].split()
                    self.wall_face[ln[0]] = float(ln[1])*(1e-3)*math.pi*self.Dia/12.0

    def facelist_gen(self):
        for psr in self.nodes.keys():
            face_string = self.nodes[psr]
            face_string2 = face_string+[face_string[0]]

            self.facelist[psr] = []
            for i in range(len(face_string)):
                self.facelist[psr].append(face_string2[i]+face_string2[i+1])

        print self.facelist

    def common_face(self):
        """
        Determines common faces between reactors and correlates the alphabetical order of the face id with the
        sign of velocity flux to determine the direction fo flow between reactors. The 'from' reactor is registered at
        index 0 and the 'to' reactor  is registered at index 1.
        :return:
        """
        for face in self.velflux.keys():
            self.commonfaces[face] = [0,0]
            for psr in self.facelist.keys():
                faces = self.facelist[psr]
                face_match = (face in faces)
                rev_face = face[::-1]
                rev_face_match = (face[::-1] in faces)
                # to order a face named as a1b2 to b2a1 to check for reversed index
                if len(face) > 2:
                    # ASCII values from 1-9 is 49-57
                    for i in range(len(face)):
                        asc_val = ord(face[i])
                        if asc_val > 48 and asc_val < 58:
                            rev_face = face[i-1] + face[i]
                            break
                    face2 = face.replace(rev_face,'')
                    if i - 1 == 0:
                        rev_face = face2+rev_face
                    else:
                        rev_face = rev_face+face2
                    rev_face_match = rev_face in faces
                    cond3 = rev_face_match and self.velflux[face] < 0.0
                    condx = cond3

                # faces defined in clockwise direction such that an outflow will have a positive value
                # for a face defined in the right order
                cond1 = face_match and self.velflux[face] > 0.0
                cond2 = face_match and self.velflux[face] < 0.0
                cond3 = rev_face_match and self.velflux[face] < 0.0
                cond4 = rev_face_match and self.velflux[face] > 0.0
                cond5 = face_match and self.velflux[face] == 0.0
                cond6 = rev_face_match and self.velflux[face] == 0.0
                cond7 = (face_match or rev_face_match) and (face in self.wall_face or rev_face in self.wall_face)
                if cond7:
                    try:
                        self.reactor_wall[psr] = [self.wall_face[face], self.reactors[psr][2]]
                    except:
                        self.reactor_wall[psr] = [self.wall_face[rev_face], self.reactors[psr][2]]
                if cond1 or cond3:
                    # outflow from current reactor
                    self.commonfaces[face][0]=psr
                elif cond2 or cond4:
                    # inflow in current reactor
                    self.commonfaces[face][1]=psr
                elif cond5 or cond6:
                    try:
                        self.zeroadjustface[face].append(psr)
                    except:
                        self.zeroadjustface[face] = [psr]
                else:

                    pass

        print self.commonfaces
        print self.zeroadjustface

    def duplicates(self,itemlist, item):
        return [i for i, x in enumerate(itemlist) if x == item]

    def psr_def(self):
        for face in self.commonfaces.keys():
            rel0 = self.commonfaces[face]
            rel= [int(x) for x in rel0]
            if (rel[0] == 0 or rel[1] == 0) and not(face in self.zeroadjustface):
                continue
            if face in self.zeroadjustface:
                # to adjust for faces with zero flux as they were outside the PIV Field of View
                rel = [int(x) for x in self.zeroadjustface[face]]
                velflux_det =0
                for f in self.facelist[rel[0]]:
                    try:
                        velflux_det += self.velflux[f]
                    except:
                        try:
                            velflux_det += self.velflux[f[::-1]]
                        except:
                            pass
                        pass
                self.velflux[face]=velflux_det
            try:
                if rel[1] in self.connect[rel[0]]:
                    ind = self.connect[rel[0]].index(rel[1])
                    flux_add = abs(self.velflux[face])
                    self.alpha_mat[rel[0]][ind] = self.alpha_mat[rel[0]][ind]+flux_add
                else:
                    self.connect[rel[0]].append(rel[1])
                    self.alpha_mat[rel[0]].append(abs(self.velflux[face]))
                #self.connect[rel[0]] = list(set(self.connect[rel[0]]))
            except:
                try:
                    self.connect[rel[0]]=[rel[1]]
                    self.alpha_mat[rel[0]]=[abs(self.velflux[face])]
                except:
                    pass

        for a in self.alpha_mat:
            vel_sum = np.sum(self.alpha_mat[a])
            #if a in self.outlet:
            #    vel_sum = np.abs(self.velflux['tu'])
            self.alpha_mat[a] = list(self.alpha_mat[a]/vel_sum)
            #if a in self.outlet:
            #    out_frac = 1-np.sum(self.alpha_mat[a])
             #   self.outlet[a] = out_frac
        for out_react in self.outlet:
            self.alpha_mat[out_react] = [1.0]
            self.outlet[out_react] = 1.0
        print self.connect
        print self.alpha_mat

    def main(self,case_num):
        case = {1:"phi06",2:"phi08",3:"phi08_N2"}
        os.chdir('C:/Users/rishikeshsampa/Documents/Research/AGNES_proletariat/'+case[case_num]+'/')
        self.readata()
        self.facelist_gen()
        self.common_face()
        self.psr_def()

if __name__=="__main__":
    obj = CRNPIV_data()
    obj.main(3)
    """obj.readata()
    obj.facelist_gen()
    obj.common_face()
    obj.psr_def()"""