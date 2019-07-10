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
Lambda = 1.0
Beta = 20

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
Bubble = 0.15916  # 2D, Beta=10, rs=1

folder = "./Beta{0}_rs{1}_lambda{2}/".format(Beta, rs, Lambda)

files = os.listdir(folder)
Num = 0
Data = None
ExtMomBin = None
AngleBin = None
ScaleBin = None
for f in files:
    if re.match("vertex"+"_pid[0-9]+.dat", f):
        print f
        with open(folder+f, "r") as file:
            line1 = file.readline()
            line2 = file.readline()
            ScaleBin = np.fromstring(line2.split(":")[1], sep=' ')
            line3 = file.readline()
            AngleBin = np.fromstring(line3.split(":")[1], sep=' ')
            line4 = file.readline()
            ExtMomBin = np.fromstring(line4.split(":")[1], sep=' ')
        Num += 1
        d = np.loadtxt(folder+f)
        if Data is None:
            Data = d
        else:
            Data += d

AngleBinSize = len(AngleBin)
ScaleBinSize = len(ScaleBin)-1
ExtMomBinSize = len(ExtMomBin)
ExtMomBin /= kF

print "Found {0} files.".format(Num)
Data /= Num

Data = Data.reshape((ScaleBinSize+1, AngleBinSize, ExtMomBinSize))

qData = np.array(Data)
# qData = np.mean(qData, axis=1)/8.0/np.pi
qData = np.mean(qData, axis=1)

verData=np.zeros(len(ExtMomBin))
for i in range(len(ExtMomBin)):
    verData[i]=integrate.simps(qData[:-1,i], ScaleBin[:-1])


def ErrorPlot(p, x, d, color, marker, label=None, size=4, shift=False):
    p.plot(x, d, marker=marker, c=color, label=label,
           lw=1, markeredgecolor="None", linestyle="--", markersize=size)
    # p.errorbar(data[:,0],data[:,1], yerr=data[:,2], c=color, ecolor=color, capsize=0, linestyle="None")
    # p.fill_between(data[:,0], data[:,1]-data[:,2], data[:,1]+data[:,2], alpha=0.5, facecolor=color, edgecolor=color)


w = 1-0.429

fig, ax = plt.subplots()
# ax=fig.add_axes()
# ax = fig.add_subplot(122)

# plt.subplot(1,2,2)
ColorList = ['k', 'r', 'b', 'g', 'm', 'c']
ColorList = ColorList*10

# ErrorPlot(ax, ScaleBin, qData[:, 20],
#             ColorList[0], 's', "Q {0}".format(ExtMomBin[20]))
# ErrorPlot(ax, ScaleBin, qData[:, 15],
#             ColorList[1], 's', "Q {0}".format(ExtMomBin[15]))

for i in range(8):
# for i in range(ScaleBinSize+1):
    # ErrorPlot(ax, ExtMomBin, qData[i, :],
    #           ColorList[i], 's', "Scale {0}".format(ScaleBin[i]))
    ErrorPlot(ax, AngleBin, Data[i, :, 0],
              ColorList[i], 's', "Order {0}".format(i))
    # ErrorPlot(ax, DataAtOrder[o], ColorList[i], 's', "Order {0}".format(o))

# ErrorPlot(ax, ExtMomBin, verData, 'y', 'o', "Sum")

# ErrorPlot(ax, Data[1][0], 'k', 's', "Diag 1")
# ErrorPlot(ax, tmp, 'm', 's', "Diag 3+c 1")
# ErrorPlot(ax, DataAll[3], 'k', 'o', "Order 3")

# ErrorPlot(ax, Data[2][1], 'g', 'o', "Order 3 counterbubble 1")
# ErrorPlot(ax, Data[2][2], 'g', '*', "Order 3 counterbubble 2")
# ErrorPlot(ax, Data[2][3], 'g', '>', "Order 3 counterbubble 3")

# ErrorPlot(ax, Data[1][1], 'olive', 'o', "Order 3 shift 1")
# ErrorPlot(ax, Data[1][2], 'olive', '*', "Order 3 shift 2")

# ErrorPlot(ax, Data[3][0], 'k', 's', "Diag 1", shift=True)
# ErrorPlot(ax, Data[3][1], 'g', 's', "Diag 2", shift=True)
# ErrorPlot(ax, Data[3][2], 'r', '*', "Diag 3", shift=True)
# ErrorPlot(ax, Data[3][3], 'b', 's', "Diag 4", shift=True)
# ErrorPlot(ax, Data[3][4], 'olive', '*', "Diag 5", shift=True)
# ErrorPlot(ax, Data[3][5], 'm', 's', "Diag 6", shift=True)
# ErrorPlot(ax, Data[3][6], 'c', '*', "Diag 7", shift=True)

# ErrorPlot(ax, Data[5], 'g', 's', "Diag 6")


x = np.arange(0, 3.0, 0.001)
y = x*0.0+Bubble
for i in range(len(x)):
    if x[i]>2.0:
        y[i]=Bubble*(1-np.sqrt(1-4/x[i]**2))

z=1.0/(1.0+y)
y=1.0-y

# y=1.0/(x*x*kF*kF+1.0)

# x = np.arange(0, 3.0, 0.001)
# y = x*0.0+Bubble
# for i in range(len(x)):
#     if x[i]>2.0:
#         y[i]=Bubble*(1-np.sqrt(1-4/x[i]**2))

# y=1.0/(x*x*kF*kF+1.0+y)

# ax.plot(x,y,'k-', lw=2)
# ax.plot(x,z,'b-', lw=2)

# ax.set_xlim([0.0, ScaleBin[-1]])
# ax.set_xlim([0.0, ExtMomBin[-1]])
ax.set_xlim([0.0, AngleBin[-1]])
# ax.set_xticks([0.0,0.04,0.08,0.12])
# ax.set_yticks([0.35,0.4,0.45,0.5])
# ax.set_ylim([-0.02, 0.125])
# ax.set_ylim([0.07, 0.125])
ax.set_xlabel("$q/k_F$", size=size)
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
