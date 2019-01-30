import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
import sys
import glob, os, re
mat.rcParams.update({'font.size': 16})
mat.rcParams["font.family"] = "Times New Roman"
size=12

rs=1.0
Lambda=1.0
Beta=25

##############   3D    ##################################
kF=(9.0*np.pi/4.0)**(1.0/3.0)/rs #3D
###### Bare Green's function    #########################
# Bubble=0.08871  # 3D, Beta=0.5, rs=1
# Bubble=0.0971916  #3D, Beta=10, rs=1
# Bubble=0.0971613  #3D, T=0.04Ef, rs=1
# Bubble= 0.097226 # 3D, zero temperature, rs=1
###### Fock dressed Green's function ###################
Bubble, Density=0.088883,0.2387 #3D, Beta=0.1, rs=1 

##############   2D    ##################################
###### Bare Green's function    #########################
# kF=np.sqrt(2.0)/rs #2D
# Bubble=0.11635  #2D, Beta=0.5, rs=1
# Bubble=0.15916  #2D, Beta=10, rs=1

ScanOrder=[1,2, 3]
# ScanOrder=[3]
Index={}
Index[1]=[1, ]
Index[2]=[1, ]
Index[3]=[1,2,3,4,5]
DataAll={}
Data={}
DataAtOrder={}
Normalization=1


folder="./Beta{0}_rs{1}_lambda{2}_eqTime/Data/".format(Beta, rs, Lambda) 
# os.chdir(folder)
files=os.listdir(folder)
for order in ScanOrder:
    Num=0
    data0=None
    for f in files:
        if re.match("Diag"+str(order)+"_[0-9]+.dat", f):
            print f
            Num+=1
            d=np.loadtxt(folder+f)
            if data0 is None:
                data0=d
            else:
                data0[:,1:]+=d[:,1:]

    print "Found {0} files.".format(Num)
    data0[:,1:]/=Num

    DataAll[order]=np.array(data0)

    Data[order]=[]
    for i in Index[order]:
        Num=0
        data=None
        for f in files:
            if re.match("Diag"+str(order)+"_[0-9]+_"+str(i)+".dat", f):
                print f
                Num+=1
                d=np.loadtxt(folder+f)
                # print f, d[0,1]
                if data is None:
                    data=d
                else:
                    data[:,1:]+=d[:,1:]
        print "Found {0} files.".format(Num)
        data[:,1:]/=Num
        Data[order].append(np.array(data))

Normalization=0.0
Num=0
for i in range(len(Data[1][0][:,0])):
    # print Data[1][0][i,0]
    if Data[1][0][i,0]>5.0*kF:
        Normalization+=Data[1][0][i,1]
        Num+=1
Normalization/=Num
Normalization/=Density

for key in DataAll.keys():
    DataAll[key][:,1]/=Normalization
    for i in range(len(Data[key])):
        Data[key][i][:,1]/=Normalization

#Additional sign due to the order
DataAll[1][:,1]*=1.0
DataAll[2][:,1]*=-1.0
DataAll[3][:,1]*=1.0

DataAtOrder[1]=np.copy(DataAll[1])
DataAtOrder[2]=np.copy(DataAll[1])
DataAtOrder[2][:,1]+=DataAll[2][:,1]
DataAtOrder[3]=np.copy(DataAll[1])
DataAtOrder[3][:,1]+=DataAll[2][:,1]
DataAtOrder[3][:,1]+=DataAll[3][:,1]


def ErrorPlot(p, d, color, marker, label=None, size=4):
    data=np.array(d)
    data[:,0]/=kF
    p.plot(data[:,0],data[:,1],marker=marker,c=color, label=label,lw=1, markeredgecolor="None", linestyle="--", markersize=size)
    # p.errorbar(data[:,0],data[:,1], yerr=data[:,2], c=color, ecolor=color, capsize=0, linestyle="None")
    # p.fill_between(data[:,0], data[:,1]-data[:,2], data[:,1]+data[:,2], alpha=0.5, facecolor=color, edgecolor=color)

w=1-0.429

fig, ax = plt.subplots()
# ax=fig.add_axes()
# ax = fig.add_subplot(122)

# plt.subplot(1,2,2)
ColorList=['k','r', 'b', 'g', 'm', 'c']

for i in range(0, len(ScanOrder)):
    o=ScanOrder[i]
    ErrorPlot(ax, DataAtOrder[o], ColorList[i], 's', "Order {0}".format(o))

# ErrorPlot(ax, DataAll[3], 'b', 'o', "Order 3")
# ErrorPlot(ax, Data[2][2], 'olive', 's', "Order 3, bubble counter 1")
# ErrorPlot(ax, Data[2][3], 'olive', 'o', "Order 3, bubble counter 2")
# ErrorPlot(ax, Data[1][1], 'm', 's', "Order 3, shift 1")
# ErrorPlot(ax, Data[1][2], 'm', 'o', "Order 3, shift 2")
# ErrorPlot(ax, Data[3][3], 'm', 's', "Diag 4")
# ErrorPlot(ax, Data[3][4], 'c', '*', "Diag 5")

# ErrorPlot(ax, Data[5], 'g', 's', "Diag 6")


x=np.arange(0,Data[1][0][-1,0]/kF,0.001)
y=x*0.0

ax.plot(x,y,'k-', lw=1)

ax.set_xlim([0.0, Data[1][0][-1,0]/kF])
# ax.set_xticks([0.0,0.04,0.08,0.12])
# ax.set_yticks([0.35,0.4,0.45,0.5])
# ax.set_ylim([0.0, 0.12])
ax.set_xlabel("$q/k_F$", size=size)
# ax.xaxis.set_label_coords(0.97, -0.01)
# # ax.yaxis.set_label_coords(0.97, -0.01)
# ax.text(-0.012,0.52, "$-I$", fontsize=size)
ax.set_ylabel("$-P(\omega=0, q)$", size=size)

# ax.text(0.02,0.47, "$\\sim {\\frac{1}{2}-}\\frac{1}{2} {\\left( \\frac{r}{L} \\right)} ^{2-s}$", fontsize=28)

plt.legend(loc=1, frameon=False, fontsize=size)
# plt.title("2D density integral")
plt.tight_layout()

# plt.savefig("spin.pdf")
plt.show()

