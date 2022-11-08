import numpy as np
import utility as ut
from const import *
from monty.serialization import loadfn

class defect:
    def __init__(self,volume,nsize,path,d_vbm=0,d_cbm=0):
        self.volume = volume
        self.nsize = nsize
        self.d_vbm = d_vbm
        self.d_cbm = d_cbm
        self.summary = self.read_yaml(path)
        

    def read_yaml(self, path):
        return loadfn(path + 'defect_energy_summary.json')

    def energy(self,efermi,temp, chem_pot_label):
        kinds = self.summary.defect_energies.keys()
        vbm = self.summary.supercell_vbm
        cbm = self.summary.supercell_cbm
        cbm = cbm - vbm
        vbm = 0
        carries_defect = 0
        rel_chem_plots =  self.summary.rel_chem_pots[chem_pot_label]

        for index,item in enumerate(kinds):
            
            each_defect = self.summary.defect_energies[item]
            each_charge = each_defect.charges
            each_energy = [ item.formation_energy for item in each_defect.defect_energies]
            each_correct = [ sum([ subitem for subitem in item.energy_corrections.values()]) \
                                    for item in each_defect.defect_energies]
            each_charge = np.array(each_charge)
            each_energy = np.array(each_energy)
            each_correct = np.array(each_correct)
            formation = each_energy + each_correct \
                                + (efermi + self.d_vbm)*each_charge 

            reservoir = sum([-diff * rel_chem_plots[elem]
                                  for elem, diff in each_defect.atom_io.items()])
            # print(reservoir)

            # Check_min_defect

            # min_index = np.argmin(formation)
            # print(str(item) + ": " + str(each_charge[min_index]))
            # print(formation[min_index]+reservoir)

            each_con = 0
            for each_charge_form,charge in zip(formation,each_charge):
                # print(each_charge_form + reservoir)
                each_con += ut.boltz_dist(self.nsize[index],self.volume,\
                    (each_charge_form + reservoir),temp)*charge
            carries_defect = carries_defect+each_con
            # print(item+': ' + '%.4e' %each_con )
            # print('Eform: ' + '%.4e' %Eform)

        return carries_defect