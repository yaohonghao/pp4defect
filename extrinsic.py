import numpy as np
import math
import utility as ut
from const import *
from monty.serialization import loadfn

class defect:
    def __init__(self,volume,nsize):
        self.volume = volume
        self.nsize = nsize
        
        self.des = 0

    def read_yaml(self, path):
        return loadfn(path + 'defect_energy_summary.json')

    def charged_defect(self,efermi,temp):
        kinds = self.des.defect_energies.keys()
        vbm = self.des.supercell_vbm
        cbm = self.des.supercell_cbm
        cbm = cbm - vbm
        vbm = 0
        carries_defect = 0

        for index,item in enumerate(kinds):
            each_defect = self.des.defect_energies[item]
            each_charge = each_defect.charges
            each_energy = [ item.formation_energy for item in each_defect.defect_energies]
            each_charge = np.array(each_charge)
            each_energy = np.array(each_energy)
            formation = each_energy + efermi*each_charge
            print(formation)
            min_index = np.argmin(formation)
            print(item)
            each_con = ut.boltz_dist(self.nsize[index],self.volume,formation[min_index],temp)*each_charge[min_index]
            carries_defect = carries_defect+each_con

        return carries_defect