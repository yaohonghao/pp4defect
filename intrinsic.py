import numpy as np


class DOS:
    def __init__(self,path,efermi,gap=0):
        self.path= path+"TDOS.dat"
        self.gap = gap

        self.dos_x, self.dos_y = self.read_dos()
        self.dos_x = self.dos_x-efermi
        self.delta = self.dos_x[1]-self.dos_x[0]
        self.VBM,self.CBM = self.find_band_edge()
        if gap == 0:
            pass
        else:
            self.scissor()

        return
         
    def read_dos(self):
        
        with open(self.path) as dos_file:
            dos=[]
            energy=[]
            while True:
                lines = dos_file.readline()
                if not lines:
                    break
                    pass
                En_temp, DOS_temp = [float(i) for i in lines.split()]
                energy.append([En_temp])
                dos.append([DOS_temp])
            energy = np.array(energy)
            dos = np.array(dos)

        return energy, dos

    def find_band_edge(self):

        start_point = np.where(self.dos_x < 0.0, self.dos_x, -np.inf).argmax() # https://www.cnpython.com/qa/170441

        while self.dos_y[start_point] != 0.0:
            start_point += 1
        VBM = start_point

        while self.dos_y[start_point] == 0.0:
            start_point += 1
        CBM = start_point

        # print("VBM: "+str(VBM)+" and the E is "+str(self.dos_x[VBM]))
        # print("CBM: "+str(CBM)+" and the E is "+str(self.dos_x[CBM]))
        return  VBM, CBM

    def scissor(self):
        """
        TODO
        shift the band edge (CBM) to more accurate value or exp. gap
        """
        return