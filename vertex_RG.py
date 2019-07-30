import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib as mat
import sys
import glob
import os
import re
mat.rcParams.update({'font.size': 16})
mat.rcParams["font.family"] = "Times New Roman"
size = 12

rs = 1.0
Lambda = 2
Beta = 20
# XType = "Tau"
XType = "Mom"
# XType = "Angle"

##############   3D    ##################################
# kF = (9.0*np.pi/4.0)**(1.0/3.0)/rs  # 3D
###### Bare Green's function    #########################
# Bubble=0.08871  # 3D, Beta=0.5, rs=1
# Bubble=0.0971916  #3D, Beta=10, rs=1
# Bubble = 0.0971613  # 3D, T=0.04Ef, rs=1
# Bubble= 0.097226 # 3D, zero temperature, rs=1
###### Fock dressed Green's function ###################
# Bubble, Density = 0.088883, 0.2387  # 3D, Beta=0.1, rs=1, Lambda=1.0

##############   2D    ##################################
###### Bare Green's function    #########################
kF = np.sqrt(2.0)/rs  # 2D
# Bubble=0.11635  #2D, Beta=0.5, rs=1
# Bubble = 0.15916/2  # 2D, Beta=10, rs=1
Bubble = 0.0795775  # 2D, Beta=20, rs=1


folder = "./Beta{0}_rs{1}_lambda{2}/".format(Beta, rs, Lambda)

files = os.listdir(folder)
Num = 0
Data0 = None
ExtMomBin = None
AngleBin = None
TauBin = None
for f in files:
    if re.match("vertex"+"_pid[0-9]+.dat", f):
        print f
        with open(folder+f, "r") as file:
            line1 = file.readline()
            line2 = file.readline()
            TauBin = np.fromstring(line2.split(":")[1], sep=' ')
            line3 = file.readline()
            AngleBin = np.fromstring(line3.split(":")[1], sep=' ')
            line4 = file.readline()
            ExtMomBin = np.fromstring(line4.split(":")[1], sep=' ')
        Num += 1
        d = np.loadtxt(folder+f)
        if Data0 is None:
            Data0 = d
        else:
            Data0 += d

AngleBinSize = len(AngleBin)
TauBinSize = len(TauBin)
ExtMomBinSize = len(ExtMomBin)
ExtMomBin /= kF

Data = None
for f in files:
    if re.match("vertex2"+"_pid[0-9]+.dat", f):
        print f
        d = np.loadtxt(folder+f)
        if Data is None:
            Data = d
        else:
            Data += d

print "Found {0} files.".format(Num)
Data /= Num
Data0 /= Num

Data = Data.reshape((AngleBinSize, ExtMomBinSize, TauBinSize))
Data0 = Data0.reshape((AngleBinSize, ExtMomBinSize, TauBinSize))

qData = np.array(Data)
qData0 = np.array(Data0)
# qData = np.sum(qData, axis=1)/AngleBinSize*2.0*np.pi
qData = np.mean(qData, axis=0)
qData0 = np.mean(qData0, axis=0)
# qData = np.mean(qData, axis=1)*2

# verData=np.zeros(len(ExtMomBin))
# for i in range(len(ExtMomBin)):
#     verData[i]=integrate.simps(qData[:,i], ScaleBin[:])

# y=np.power(ScaleBin[:-1], 3)
# print ScaleBin[-5:-1]
# print ScaleBin
# print integrate.simps(y, ScaleBin[:-1])


def ErrorPlot(p, x, d, color, marker, label=None, size=4, shift=False):
    p.plot(x, d, marker=marker, c=color, label=label,
           lw=1, markeredgecolor="None", linestyle="--", markersize=size)


w = 1-0.429

fig, ax = plt.subplots()
# ax=fig.add_axes()
# ax = fig.add_subplot(122)

# plt.subplot(1,2,2)
ColorList = ['k', 'r', 'b', 'g', 'm', 'c', 'navy', 'y']
ColorList = ColorList*40

if(XType == "Scale"):
    for i in range(ExtMomBinSize/4):
        index = 4*i
        ErrorPlot(ax, ScaleBin[:-2], diffData[:-2, index],
                  ColorList[i], 's', "Q {0}".format(ExtMomBin[index]))
    ax.set_xlim([0.0, ScaleBin[-2]])
    ax.set_xlabel("$Scale$", size=size)
elif(XType == "Tau"):
    for i in range(ExtMomBinSize/8):
        index = 8*i
        ErrorPlot(ax, TauBin, qData0[index, :],
                  ColorList[2*i], 's', "Order 1, Q {0}".format(ExtMomBin[index]))
        ErrorPlot(ax, TauBin, qData[index, :],
                  ColorList[2*i+1], 's', "Order 2, Q {0}".format(ExtMomBin[index]))
    ax.set_xlim([0.0, TauBin[-1]+0.1])
    ax.set_xlabel("$Tau$", size=size)
elif (XType == "Mom"):

    # qData = np.sum(qData, axis=1)/(128.0/np.pi)
    # qData0 = np.sum(qData0, axis=1)*Beta/kF**2/TauBinSize/40.65
    qData = np.sum(qData, axis=1)*Beta/kF**2/TauBinSize
    qData0 = np.sum(qData0, axis=1)*Beta/kF**2/TauBinSize
    # qData0 = qData0/qData0[0]*12.50
    qData0 = 8.0*np.pi/(ExtMomBin**2*kF**2+Lambda)-qData0
    # qData=8.0*np.pi/(ExtMomBin**2*kF**2+Lambda)-qData

    ErrorPlot(ax, ExtMomBin, qData0[:],
              ColorList[1], 's', "Order 2")

    ErrorPlot(ax, ExtMomBin, qData[:],
              ColorList[2], 's', "Order 3")
    # for i in range(ScaleBinSize/32):
    #     # print i, index
    #     # print ScaleBin[index]
    #     index = 32*i
    # # for i in range(ScaleBinSize+1):
    #     ErrorPlot(ax, ExtMomBin, qData[index, :],
    #               ColorList[i], 's', "Scale {0}".format(ScaleBin[index]))
    # ErrorPlot(ax, AngleBin, Data[i, :, 8],
    #           ColorList[i], 's', "Order {0}".format(i))
    # ErrorPlot(ax, DataAtOrder[o], ColorList[i], 's', "Order {0}".format(o))

    # ErrorPlot(ax, ExtMomBin, verData, 'y', 'o', "Sum")

    # x = np.arange(0, 3.0, 0.001)
    # y = x*0.0+Bubble
    # for i in range(len(x)):
    # if x[i]>2.0:
    # y[i]=Bubble*(1-np.sqrt(1-4/x[i]**2))

    # z=1.0/(1.0+y)
    # y=1.0-y

    # y=1.0/(x*x*kF*kF+1.0)

    x = np.arange(0, 3.0, 0.001)
    y = x*0.0+Bubble
    for i in range(len(x)):
        if x[i] > 2.0:
            y[i] = Bubble*(1-np.sqrt(1-4/x[i]**2))
    y0 = 8.0*np.pi/(x*x*kF*kF+Lambda)
    # ym=y0-y0*y0*y
    yphy = 8.0*np.pi/(x*x*kF*kF+Lambda+y*8.0*np.pi)

    # ax.plot(x, yphy, 'k-', lw=2, label="physical")
    ax.plot(x, y0, 'k-', lw=2, label="original")

    # ax.plot(x, y0*y0*y, 'r-', lw=2, label="wrong")

    ax.set_xlim([0.0, ExtMomBin[-1]])
    ax.set_xlabel("$q/k_F$", size=size)

elif(XType == "Angle"):
    for i in range(ScaleBinSize/8):
        # print i, index
        # print ScaleBin[index]
        index = 8*i
        ErrorPlot(ax, AngleBin, Data[index, :, 5],
                  ColorList[i], 's', "Q {0}".format(ScaleBin[index]))
    ax.set_xlim([0.0, AngleBin[-1]])
    ax.set_xlabel("$Angle$", size=size)
# ax.set_xticks([0.0,0.04,0.08,0.12])
# ax.set_yticks([0.35,0.4,0.45,0.5])
# ax.set_ylim([-0.02, 0.125])
# ax.set_ylim([0.07, 0.125])
# ax.xaxis.set_label_coords(0.97, -0.01)
# # ax.yaxis.set_label_coords(0.97, -0.01)
# ax.text(-0.012,0.52, "$-I$", fontsize=size)
ax.set_ylabel("$-\Gamma_4(\omega=0, q)$", size=size)

# ax.text(0.02,0.47, "$\\sim {\\frac{1}{2}-}\\frac{1}{2} {\\left( \\frac{r}{L} \\right)} ^{2-s}$", fontsize=28)

plt.legend(loc=1, frameon=False, fontsize=size)
# plt.title("2D density integral")
plt.tight_layout()

# plt.savefig("spin_rs1_lambda1.pdf")
plt.show()
