import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import fft, ifft
import cmath
#from scipy.interpolate import spline
#file = sys.argv[1]

###############################################################################################
#
# Read velocities from trajectory
#
###############################################################################################
atom_num=384   #total atom number
frame_num=20000  #total frame number and also the FFT length
t_total=10000  #total time in fs
T=300
amu=1.660538921e-27
kb=1.38064852e-23
ang=1e-10
fs=1e-15
angbyfs=0.529177/0.024189   #velocity unit conversion from bohr/aut to Ang/fs

traj=sys.argv[1]
frame_segment=np.zeros((atom_num,3))
g=np.zeros((atom_num,frame_num*3))
count=l=0
with open (traj, 'rt') as input1:  #read the velocity trajectory
    for line in input1:
        l=l+1
        if l>atom_num+2:
            l=1
            count=count+1
        if l<3:
            continue
        wordList=[]
        for word in line.split():
            wordList.append(word)
        g[l-3,count*3+0]=float(wordList[1])
        g[l-3,count*3+1]=float(wordList[2])
        g[l-3,count*3+2]=float(wordList[3])
        if count==frame_num-1:
            break
mass=np.zeros((atom_num*3,1))

i=0
with open('mass','rt') as input2:
    for line in input2:
        wordList=[]
        for word in line.split():
            wordList.append(word)
        mass[i,0]=float(wordList[0])       
        mass[i+1,0]=float(wordList[0])
        mass[i+2,0]=float(wordList[0])
        i=i+3

h=np.transpose(g)
single_v=np.zeros((frame_num, 3))
total_v=np.zeros((frame_num,0))
sv_total=0

for i in range(atom_num):
    for j in range(frame_num*3):
        single_v[int(j/3),j%3]=h[j,i]
    total_v=np.append(total_v,single_v,1)  #this is a frame_num*3atom_num matrix

total_v=total_v*angbyfs

v_filter=np.zeros((frame_num,atom_num*3),dtype=complex)
sv_total=np.zeros((0, frame_num))
fft_factor=t_total/(frame_num)**2    #FFT length transformation 
for i in range(atom_num*3):
    sv_single=fft(total_v[:,i])
    v_filter[:,i]=sv_single
    for j in range(frame_num):
        sv_single[j]=abs(sv_single[j])
        sv_single[j]=(sv_single[j])**2*fft_factor
        sv_single[j]=sv_single[j]*mass[i,0]
    sv_total=np.append(sv_total, np.array([sv_single]), 0)   #this is a 3atom_num*frame_num matrix

dos=np.zeros((frame_num,3))
dens=open('DoS', 'w')

conversion=ang**2/fs**2*amu/kb
twobykbT=2/T
for i in range(frame_num):
    dos[i,0]=i
    y=np.sum(sv_total[:,i])
    dos[i,1]=y
    dos[i,1]=dos[i,1]*conversion*twobykbT
    dens.write(str(dos[i,0])+'   '+str(dos[i,1])+'\n')
