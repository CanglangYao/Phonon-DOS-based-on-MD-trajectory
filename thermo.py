import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import fft, ifft
h= 6.62607015e-19  #planck constant J*fs
kb= 1.38064852e-23
T=300
j2eV=6.242e18
atom_num=384
frame_num=20000
t_total=10000.0  #total time in fs
angbyfs=0.529177/0.024189  #velocity unit in this code is A/fs, we'll return to the cp2k unit
i=0
mode=1 #input the vibration mode number 
g=np.zeros((int(frame_num), 3))
with open ('DoS', 'rt') as f1:
    for line in f1:
        wordList = []
        for word in line.split():
            wordList.append(word)
        g[i,0]=(i+1)/t_total
        g[i,1]=float(wordList[1])
        bhv=h*g[i,0]/(kb*T)
        ws=bhv/(math.exp(bhv)-1)-math.log(1-math.exp(-bhv))
        g[i,2]=g[i,1]*ws*1/t_total
        i=i+1
        if i==int(frame_num/2):       #FFT also generates imaginary result, which has same strength with thre real one, we discard the imaginary value.
            break
entropy=np.sum(g[:,2])*kb
e_vib=-T*entropy*j2eV
e_vib=e_vib/32
print('vibrational energy',e_vib)

