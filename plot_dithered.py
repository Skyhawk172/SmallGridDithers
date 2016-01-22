import numpy as np
import matplotlib.pyplot as plt
import glob, os, sys



dir=sys.argv[1]
ditherstep=dir[-2:]
os.chdir('/Users/lajoie/Documents/Work/Projects/JWST/Simulations/Coronagraphs/MIRI/Dither-LOCI/Results/'+dir)
#os.chdir('./'+dir)

files=glob.glob('run*.txt')
print files
nx=1
ny=5

ax=[]
f, ax = plt.subplots(nx,ny, sharex=True, sharey=True, figsize=(12,8))
f.subplots_adjust(wspace=0)
f.subplots_adjust(hspace=0)

fontsize=22
markersize=15

nf=0
for i in xrange(nx):
    for j in xrange(ny):   
        pos =np.loadtxt(files[nf],unpack=True)
        pos2=np.reshape(pos,(11,2))
        pos3=np.delete(pos2,1,0)   

        if nx==1: 
            ax[j].plot(pos3[:,0],pos3[:,1],'.' ,markersize=markersize)
            ax[j].plot(pos3[0,0],pos3[0,1],'r.',markersize=markersize)
            ax[j].set_xlim(-50,50)
            ax[j].set_ylim(-50,50)
            ax[j].axhline(color='k',linestyle="--")
            ax[j].axvline(color='k',linestyle="--")
            ax[j].text(10., 40., files[nf][:-4],fontsize=fontsize)
            if j==0: ax[j].set_ylabel('y (mas)',fontsize=fontsize)
            ax[j].set_xlabel('x (mas)',fontsize=fontsize)
            if i==0 and j==5: ax[j].set_title(sys.argv[1])
            ax[j].yaxis.set_tick_params(labelsize=fontsize)
            ax[j].xaxis.set_tick_params(labelsize=fontsize)
            ax[j].set_aspect("equal")
  
        else:
  
            ax[i,j].plot(pos3[:,0],pos3[:,1],'.' ,markersize=markersize)
            ax[i,j].plot(pos3[0,0],pos3[0,1],'r.',markersize=markersize)
            ax[i,j].set_xlim(-50,50)
            ax[i,j].set_ylim(-50,50)
            ax[i,j].axhline(color='k',linestyle="--")
            ax[i,j].axvline(color='k',linestyle="--")

            if j==0: ax[i,j].set_ylabel('y (mas)',fontsize=fontsize)
            if j==0: ax[i,j].set_xlabel('x (mas)',fontsize=fontsize)
            ax[i,j].text(10., 40., files[nf][:-4],fontsize=fontsize)
            if i==0 and j==5: ax[i,j].set_title(sys.argv[1])
            ax[j].yaxis.set_tick_params(labelsize=fontsize)
            ax[j].xaxis.set_tick_params(labelsize=fontsize)

        nf+=1

ax[np.ceil(ny/2)].set_title("Dither step=%s mas" %ditherstep,fontsize=fontsize)

plt.savefig("/Users/lajoie/Documents/Work/Projects/JWST/CWG/SGD/try.pdf")
raw_input()



