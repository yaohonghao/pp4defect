from cProfile import label
from intrinsic import DOS
from extrinsic import defect
import utility as ut
import numpy as np
import matplotlib.pyplot as plt
from monty.serialization import loadfn
np.set_printoptions(suppress=False)

if __name__ == "__main__":
    
    """
    Pass
    Case: Fermi level in Silicon moves during increasing temperature
    # """
    # path = 'C:\\Users\\yaoho\\Desktop\\defect-SrCuSb\\'
    # dos = DOS(path=path,efermi=5.3631)
    # # BaCuSb : 4.8124
    # # BaAgSb : 4.9657
    # # SrCuSb : 5.3631
    # # SrAgSb : 5.5693
    # delta_T = 50
    # criteria = 2 # less than 10^criteria cm^-3
    # volume = 40.89 # in Å^3
    # for i in range(21):
    #     temperature = delta_T*i
    #     ut.find_fermi_intrinsic(dos,temperature,criteria,volume)
    
    
    """
    Step 1: 
    """
    # path = 'C:\\Users\\yaoho\\Desktop\\defect-BaAgSb\\' 
    # pydefect_summary = loadfn(path + 'defect_energy_summary.json')
    # for item in pydefect_summary.defect_energies.keys():
    #     print(item)
    # print("Defects Number: " + str(len(pydefect_summary.defect_energies.keys())))
    """
    Step 2: 
    """
    # SrCuSb
    path = 'C:\\Users\\yaoho\\Desktop\\defect-SrCuSb\\' 
    volume = 156.7415
    d_cbm=-0.01375
    d_vbm= 0.17425
    nsize = [4 for _ in range(6)]
    nsize.extend([2, 2, 2])
    efermi = 5.3631
    dos = DOS(path=path,efermi= efermi)
    
    # BaCuSb
    # path = 'C:\\Users\\yaoho\\Desktop\\defect-BaCuSb\\' 
    # volume = 176.27
    # d_cbm=-0.013266667
    # d_vbm= 0.164033333
    # nsize = [4 for _ in range(6)]
    # nsize.extend([2, 2, 2])
    # efermi = 4.8124
    # dos = DOS(path=path,efermi= efermi)

    # path = 'C:\\Users\\yaoho\\Desktop\\defect-BaAgSb\\' 
    # volume = 191.04
    # d_cbm=-0.021666667
    # d_vbm= 0.153433333
    # nsize = [4 for _ in range(6)]
    # nsize.extend([2, 2, 2])
    # dos = DOS(path=path,efermi= 4.9657)

    # path = 'C:\\Users\\yaoho\\Desktop\\defect-SrAgSb\\' 
    # volume = 170.92
    # d_cbm= -0.01375
    # d_vbm= 0.17425
    # nsize = [4 for _ in range(6)]
    # nsize.extend([2, 2, 2])
    # efermi= 5.5693
    # dos = DOS(path=path,efermi= efermi)
    
    delta_T = 20
    criteria = 8
    charge_defect = defect(volume,nsize,path,d_vbm=d_vbm, d_cbm=d_cbm)
    ## print(ut.find_defect_fermi(dos,charge_defect,'A',300,criteria,volume))
    
    for i in range(41):
        temperature = delta_T*i
        atomi_point = 'B'
        E_fermi=ut.find_defect_fermi(dos,charge_defect,atomi_point,\
            temperature,criteria,volume)
        con = -charge_defect.energy(E_fermi,temperature,atomi_point)
        int_con = ut.intrinsic(dos,fermi=E_fermi,temp=temperature,volume=volume)
        # print("%e"%str((con+int_con)/2))
        # print("%e"%con+' E_femri: %e'%E_fermi)
        # print(str(temperature)+ ' E_femri: %e'%E_fermi)
        # print("%e"%con)
        # print('%e'%int_con)
        print("%e"%((con+int_con)/2.))

            
    # E_fermi = 4.918750e-03
    # print( '%e'%float(-charge_defect.energy( E_fermi,300,'B')))
    # print('%e'%ut.intrinsic(dos,fermi=E_fermi,temp=300.0,volume=volume))

    # for i in range(6):
    #     ut.intrinsic(dos,fermi=-0.1*i,temp=300,volume=volume)
    #     print('%e'%ut.intrinsic(dos,fermi=-0.1*i,temp=300,volume=volume))   
 
    
    # con = charge_defect.energy(0.0,0,'D')
    # print("%e"%con)

    # int_def = defect(volume,nsize)
    # int_def.des = int_def.read_yaml(path)

    # for i in range(50):
    #     con = -charge_defect.energy(0.5,i*20,'D')
        # print("%e"% con)



    """
    TEST: Density of States Pass!
    """
    # path = 'C:\\Users\\yaoho\\Desktop\\defect-BaAgSb\\' 
    # dos = DOS(path=path,efermi= 4.9657)
    # plt.plot(dos.dos_x,dos.dos_y,label="BaAgSb")
    # path = 'C:\\Users\\yaoho\\Desktop\\defect-SrAgSb\\' 
    # dos = DOS(path=path,efermi= 5.5693)
    # plt.plot(dos.dos_x,dos.dos_y,label="SrAgSb")
    # path = 'C:\\Users\\yaoho\\Desktop\\defect-BaCuSb\\' 
    # dos = DOS(path=path,efermi= 4.8124)
    # plt.plot(dos.dos_x,dos.dos_y,label="BaCuSb")
    # path = 'C:\\Users\\yaoho\\Desktop\\defect-SrCuSb\\' 
    # dos = DOS(path=path,efermi= 5.3631)
    # plt.plot(dos.dos_x,dos.dos_y,label="SrCuSb") 
    # plt.xlim([-1,1.2])
    # plt.ylim([0,20])
    # plt.legend()
    # plt.show()
    

    """
    TEST: Fermi_dist Pass!
    """
    # fermi_range =  [0.01*i-1. for i in range(201)]
    # fermi_range = np.array(fermi_range)
    # distr_range = np.zeros_like([fermi_range for _ in range(3)])
    # for index_T,temp in enumerate([0, 300, 500]):
    #     for index_C,con in enumerate(fermi_range):
    #         print(index_C)
    #         distr_range[index_T][index_C] = ut.fermi_dist(0,con,temp)
    # plt.plot(fermi_range,distr_range[0],label="0 K")
    # plt.plot(fermi_range,distr_range[1],label="300 K")
    # plt.plot(fermi_range,distr_range[2],label="500 K")
    # plt.legend()
    # plt.show()