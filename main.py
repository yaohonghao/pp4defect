from intrinsic import DOS
import utility as ut

if __name__ == "__main__":
    
    """
    Pass
    Case: Fermi level in Silicon moves during increasing temperature
    """
    path = 'C:\\Users\\yaoho\\Desktop\\'
    dos = DOS(path=path,efermi=5.6320)
    
    delta_T = 100
    criteria = 10 # less than 10^criteria cm^-3
    volume = 40.89
    for i in range(10):
        temperature = delta_T*i
        print(ut.find_defect_fermi(dos,temperature,criteria,volume))
     # print(E_fermi)
    