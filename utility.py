import enum
import numpy as np
import math
from const import *


def fermi_int(dos,Efermi,temp):
    # con_elec = 
    start = dos.CBM
    while dos.dos_x[start] < dos.dos_x[dos.CBM]+2:
        start += 1
    upper = start
    # print(type(dos.dos_y[dos.CBM:upper]))
    # print(dos.dos_x[upper])
    num_ele = [ float(fermi_dist(Efermi,i,temp)) for i in dos.dos_x[dos.CBM:upper]]
    num_ele = np.array(num_ele)
    # print(type(num_ele))
    abc = dos.dos_y[dos.CBM:upper]
    abc = [float(item) for item in abc]
    abc = np.array(abc)
    c_elec = np.dot(num_ele,abc)#/dos.delta
    
    count = 0
    for index, ele in enumerate(dos.dos_x[dos.CBM:upper]):
        count += fermi_dist(Efermi,ele,temp) * dos.dos_y[index]
    # print(float(count * dos.delta))

    start = dos.VBM
    while dos.dos_x[start] > dos.dos_x[dos.VBM]-2:
        start -= 1
    lower = start
    # print(dos.dos_x[lower])
    num_hole = [float(1-fermi_dist(Efermi,i,temp)) for i in dos.dos_x[lower:dos.VBM]]
    num_hole = np.array(num_hole)
    # print(type(num_ele))
    abc = dos.dos_y[lower:dos.VBM]
    abc = [float(item) for item in abc]
    abc = np.array(abc)
    # print(num_ele.shape)
    # print(abc.shape)
    c_hole = np.dot(abc,num_hole)#/dos.delta

    # count_1 = 0
    # for index, ele in enumerate(dos.dos_x[lower:dos.VBM]):
    #     count_1 += (1-fermi_dist(Efermi,ele,temp)) * dos.dos_y[index]
    # print(c_hole)
    # print(c_elec)
    return -c_elec*dos.delta, c_hole*dos.delta

def boltz_dist(nsize, volume, eform, temp):
    # print(eform)
    temp = temp + 1e-3
    if eform < 0.0:
        # concentration = (nsize / volume) * math.exp(0) * (-eform/ (k_b*temp/eV))*10**21
        # print("%e"%concentration)
        
        # b =( 1 - eform/ (k_b*temp/eV))*10.**24
        a = nsize / volume
        return a*10.**24
        # return (nsize / volume) * math.exp(0) * ((-eform/ (k_b*temp/eV))*10**21)
    elif eform == 0.0:
        return (nsize / volume) *10.**24
    else:
        # print( (nsize / volume) * 10**21)
        return (nsize / volume) * math.exp(-eform/ (k_b*temp/eV))*10**24

def fermi_dist(efermi, epsilon,temp):
    temp += 0.1
    power = (epsilon - efermi)/k_b/temp*eV
    if power > 100:
        return 0
    else:
        return 1./(math.e**power+1)

def intrinsic(dos,fermi,temp,volume):
    n, p = fermi_int(dos,fermi,temp)
    # print("n-type: %e" %n)
    # print("p-type: %e" %p)
    num  = float(n+p)
    con  = compute_carrier(num,volume)
    # print("concentration: %e"%con)
    return con


def compute_carrier(num, volume):
    return float(num*10**24/volume)

def extract_power(x): 
    return np.sign(x)*np.log10(abs(x)+1)   # why adding 1 in np.log()?

def self_consist(dos, fermi_0,y_0,fermi_1,y_1,temp,volume,criteria=10):
    global count 
    fermi_2 = (fermi_1 + fermi_0) / 2
    con_2 = intrinsic(dos, fermi=fermi_2,temp=temp,volume=volume)
    y_2 = extract_power(con_2)
    # print("The fermi is:" + str(fermi_2) + " The log is:" + str(y_2))

    if abs(y_2) < criteria:
        return  fermi_2
    else:
        count += 1
        if y_0 * y_2 < 0:
            return self_consist(dos, fermi_0,y_0,fermi_2,y_2,temp=temp,volume=volume,criteria=criteria)
        else:
            return self_consist(dos, fermi_1,y_1,fermi_2,y_2,temp=temp,volume=volume,criteria=criteria)

def self_consist2(dos, Defecet,miu_atom, fermi_0,y_0,fermi_1,y_1,temp,volume,criteria=10):
    global count 
    fermi_2 = (fermi_1 + fermi_0) / 2
    con_2 = intrinsic(dos, fermi=fermi_2,temp=temp,volume=volume)
    con_2 += Defecet.energy(fermi_2,temp,miu_atom)
    y_2 = extract_power(con_2)
    # print('fermi_2: ' + str(fermi_2))
    # print("The fermi is:" + str(fermi_2) + " The log is:" + str(y_2))
    if abs(y_2) < criteria:
        return  fermi_2
    else:
        count += 1
        if y_0 * y_2 < 0:
            if abs(fermi_0-fermi_2) <= 0.005:
                return (fermi_0 + fermi_2)/2
            else:
                # print("Fail1")
                return self_consist2(dos, Defecet,miu_atom, fermi_0,y_0,fermi_2,y_2,temp=temp,volume=volume,criteria=criteria)
        else:
            if abs(fermi_1-fermi_2) <= 0.01:
                # print("Attention, this values is a compirsed result")
                return (fermi_1 + fermi_2)/2
            else:
                # print("Fail2")
                return self_consist2(dos, Defecet,miu_atom, fermi_1,y_1,fermi_2,y_2,temp=temp,volume=volume,criteria=criteria)

def find_fermi_intrinsic(dos, temp, criteria, volume):
    fermi_0 = int((dos.CBM + dos.VBM)/2)
    fermi_0 = float(dos.dos_x[fermi_0])
    con_0 = intrinsic(dos,fermi=fermi_0,temp=temp,volume=volume)
    y_0 = extract_power(con_0)
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    step = np.sign(con_0)*(dos.CBM-dos.VBM)*0.1 
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # print('==================')
    # print(str(step))
    # print('==================')
    fermi_1 = fermi_0 +  step
    con_1 = intrinsic(dos,fermi=fermi_1,temp=temp,volume=volume)
    y_1 = extract_power(con_1)
    if y_0 * y_1 <= 0:
        check = True
    else:
        check = False

    assert check ," Warning: the second point has the same sign with the first!"

    if abs(y_0) < criteria:
        E_fermi = fermi_0
    elif abs(y_1) < criteria:
        E_fermi = fermi_1
    else:
        global count 
        count = 1
        E_fermi = self_consist(dos,fermi_0,y_0,fermi_1,y_1,temp=temp,volume=volume,criteria=criteria)
        # print("During " +str(count) + " iterations")
    return E_fermi

def find_defect_fermi(dos, Defecet,miu_atom, temp, criteria, volume):
    fermi_0 = int((dos.CBM + dos.VBM)/2)
    fermi_0 = float(dos.dos_x[fermi_0])
    con_0 = intrinsic(dos,fermi=fermi_0,temp=temp,volume=volume)
    con_0 += Defecet.energy(fermi_0,temp,miu_atom)
    y_0 = extract_power(con_0)
    step = -1 #np.sign(con_0)*(dos.CBM-dos.VBM)*0.1
    fermi_1 = fermi_0 +  step
    # print('==================')
    # print('Gap: ' + str(dos.CBM-dos.VBM))
    # print('step: ' + str(step))
    # print('fermi_1: ' + str(fermi_1))
    # print('fermi_0: ' + str(fermi_0))
    # print('==================')
    
    con_1 = intrinsic(dos,fermi=fermi_1,temp=temp,volume=volume)
    con_1 += Defecet.energy(fermi_1,temp,miu_atom)
    y_1 = extract_power(con_1)

    if y_0 * y_1 <= 0:
        check = True
    else:
        check = False

    assert check ," Warning: the second point has the same sign with the first!"

    if abs(y_0) < criteria:
        E_fermi = fermi_0
    elif abs(y_1) < criteria:
        E_fermi = fermi_1
    else:
        global count 
        count = 1
        E_fermi = self_consist2(dos, Defecet,miu_atom,fermi_0,y_0,fermi_1,y_1,temp=temp,volume=volume,criteria=criteria)
        # print("During " +str(count) + " iterations")
    return E_fermi