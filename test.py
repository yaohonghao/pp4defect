from intrinsic import DOS
from extrinsic import defect
import utility as ut
import numpy as np
np.set_printoptions(suppress=False)

if __name__ == "__main__":
    
    """
    Pass
    Case: Fermi level in Silicon moves during increasing temperature
    # """
    # path = 'C:\\Users\\yaoho\\Desktop\\'
    # dos = DOS(path=path,efermi=5.6320)
    
    # delta_T = 100
    # criteria = 10 # less than 10^criteria cm^-3
    # volume = 40.89 # in Å^3
    # for i in range(10):
    #     temperature = delta_T*i
    #     print(ut.find_fermi_intrinsic(dos,temperature,criteria,volume))
     # print(E_fermi)


    
    # path = 'C:\\Users\\yaoho\\Desktop\\defect-BaCuSb\\' 
    # pydefect_summary = loadfn(path + 'defect_energy_summary.json')
    # for item in pydefect_summary.defect_energies.keys():
    #     print(item)
    nsize = [2, 2, 2, 2, 2, 2, 2, 2, 2]
    volume= 156.74
    path = 'C:\\Users\\yaoho\\Desktop\\defect-BaCuSb\\' 

    charge_defect = defect(volume,nsize,path)
    print( '%e'%float(charge_defect.energy(0.3,300,'C')))
    # int_def = defect(volume,nsize)
    # int_def.des = int_def.read_yaml(path)
    # con = int_def.charged_defect(0,300)
    # print("%e"% con)
    