#!/usr/bin/env python
import numpy as np
import unionfind
import random
from nullspace import rank, nullspace
from numpy.linalg import matrix_rank
import itertools 
import sys,os

random.seed(2)
OrderedTime={}
DUPLICATE=1000
for o in range(1,10):
    OrderedTime[o]=[]
    for dup in range(DUPLICATE):
        # timelist=range(o)
        # random.shuffle(timelist)
        timelist=[random.random() for _ in range(o-1)]
        dtimelist=[0.0,random.random()]
        for e in timelist:
            dtimelist.append(e)
            dtimelist.append(e)
        OrderedTime[o].append(dtimelist)

def ReNameDiag(DiagListInput, Initial=1):
    DiagList=list(DiagListInput)
    for d in DiagList:
        for i in range(len(d)/2):
            d[i*2]=d[i*2]*2-2+Initial*2
            d[i*2+1]=d[i*2+1]*2+1-2+Initial*2
    return DiagList

def Mirror(Index):
    if Index%2==0:
        return Index+1
    else:
        return Index-1

# def AttachExternalVertex(Diag):
    # PolarDict=dict()
    # for i in range(len(Diag)):
        # a=list(Diag)
        # a.append(Diag[i])
        # a[i]=0
        # for j in range(len(a)):
            # b=list(a)
            # b.append(a[j])
            # b[j]=1
            # c=list(b)
            # c[0]=b[-2]
            # c[1]=b[-1]
            # c[2:]=b[:-2]
            # if j<len(Diag):
                # PolarDict[tuple(c)]=(i+2,j+2)
            # else:
                # PolarDict[tuple(c)]=(i+2,0)
    # return PolarDict

def AttachExternalVertex(Diag, ZMomentum, Sym):
    PolarDict=dict()
    SHIFT=2
    order=len(Diag)/2+1
    for i in range(2, len(Diag)+2):
        #Initialization
        #d[i]<== 1 <== 0 <== i
        d=[0,1]+list(Diag)
        d[1]=d[i]
        d[0]=1
        d[i]=0

        Momentum=np.zeros([order+1, 2*order], dtype=int)
        Momentum[1:,2:]=ZMomentum
        Momentum[1:,0]=ZMomentum[:,i-SHIFT]
        Momentum[1:,1]=ZMomentum[:,i-SHIFT]
        Momentum[0,0]=1

        if not CheckConservation(d, Momentum, IsPolarization=True):
            print "Momentum does not conserve or rank is wrong!"
            sys.exit(0)

        # print "Start with: ", d
        PolarDict[tuple(d)]=[Momentum, Sym, TimeSign(d)]
        ToVisit=[d[1], Mirror(d[1])]
        StartPermu=[tuple(d), tuple(d)]
        StartMom=[Momentum, Momentum]
        Visited=[0]
        # print StartMom[0]
        while len(ToVisit)!=0:
            Index=ToVisit.pop()
            Permutation=list(StartPermu.pop())
            Mom=np.copy(StartMom.pop())
            if Index in Visited:
                continue

            if Permutation[1]!=Index and Permutation[1]!=Mirror(Index):
                print "wrong!", Permutation, Index
                sys.exit()
            #NextVertex<==1<===PreVertex, Target<======Index
            #NextVertex<=======PreVertex, Target<==1<==Index, 
            Target=Permutation[Index]
            NextVertex=Permutation[1]
            PrevVertex=Permutation.index(1)
            Permutation[1]=Target
            Permutation[PrevVertex]=NextVertex
            Permutation[Index]=1

            deltaMom=np.copy(Mom[:,PrevVertex]-Mom[:,1])
            Mom[:,1]=Mom[:,Index]
            Mom[:,Index]+=deltaMom

            if not CheckConservation(Permutation, Mom, IsPolarization=True):
                print "Momentum does not conserve or rank is wrong!"
                sys.exit(0)

            PolarDict[tuple(Permutation)]=[Mom, Sym, TimeSign(Permutation)]

            Visited.append(Index)

            if Target not in Visited:
                ToVisit.append(Target)
                ToVisit.append(Mirror(Target))
                StartPermu.append(tuple(Permutation))
                StartPermu.append(tuple(Permutation))
                StartMom.append(Mom)
                StartMom.append(Mom)
        # print len(Visited)
    return PolarDict

def TimeSign(permutation, FermiSign=1):
    order=len(permutation)/2
    timesign=[]
    for time in OrderedTime[order]:
        sign=FermiSign
        for i in range(len(permutation)):
            sign*=np.sign(time[permutation[i]]-time[i])
        timesign.append(sign)
    return np.array(timesign)

def swap(array, i, j):
    array = list(array)
    array[i], array[j] = array[j], array[i]
    return tuple(array)

def IsConnected(permutation, reference, InteractionPairs):
    diagram=set(InteractionPairs)
    for i in range(len(permutation)):
        diagram.add((reference[i], permutation[i]))

    n_node = len(InteractionPairs)*2
    diagram_union = unionfind.UnionFind(n_node)

    for edge in diagram:
        if edge[0]!=edge[1] and not diagram_union.is_connected(edge[0], edge[1]):
            diagram_union.union(edge[0], edge[1])
    return diagram_union.get_n_circles() == 1

def GetInteractionPairs(Order):
    return [(2*i,2*i+1) for i in range(Order)]

def GetReference(Order):
    return range(2*Order)

def HasTadpole(permutation, reference):
    for i in range(len(permutation)):
        if reference[i]==permutation[i]:
            return True
    return False

def HasFock(permutation, reference):
    for i in range(len(reference)):
        # end=reference[i]
        end=permutation[i]
        if i==0 or i==1:
            continue
        if abs(i-end)==1 and min(i, end)%2==0:
            return True
    return False

def swap(array, i, j):
    array = list(array)
    array[i], array[j] = array[j], array[i]
    return tuple(array)

def swap_interaction(permutation, m, n, k, l):
    permutation = list(permutation)
    mp,np,kp,lp=(permutation.index(e) for e in (m,n,k,l))
    permutation[mp]=k
    permutation[kp]=m
    permutation[np]=l
    permutation[lp]=n
    permutation[m],permutation[k]=permutation[k],permutation[m]
    permutation[n],permutation[l]=permutation[l],permutation[n]
    return tuple(permutation)

def swap_LR(permutation, i, j):
    # print permutation, i, j
    permutation = list(permutation)
    ip,jp=permutation.index(i),permutation.index(j)
    permutation[ip]=j
    permutation[jp]=i
    permutation[i],permutation[j]=permutation[j],permutation[i]
    # print "after", permutation
    return tuple(permutation)

def swap_LR_Hugen(permutation, i, j):
    permutation = list(permutation)
    permutation[i],permutation[j]=permutation[j],permutation[i]
    return swap_LR(permutation, i, j)
    # return tuple(permutation)

def swap_LR_Hugen_Backward(permutation, i, j):
    permutation = list(permutation)
    permutation[i],permutation[j]=permutation[j],permutation[i]
    return permutation
    # return tuple(permutation)

def check_Unique_Permutation(permutation, PermutationDict, TimeRotation):
    Order = len(permutation)/2
    Deformation = [permutation]

    if TimeRotation:
        for idx in range(1, Order):
            for i in range(len(Deformation)):
                for j in range(1, idx):
                    Deformation.append(swap_interaction(Deformation[i], idx*2, idx*2+1, j*2, j*2+1))

    for idx in range(1,Order):
        for i in range(len(Deformation)):
            Deformation.append(swap_LR(Deformation[i], idx*2, idx*2+1))

    for idx in range(1,Order):
        for i in range(len(Deformation)):
            Deformation.append(swap_LR_Hugen(Deformation[i], idx*2, idx*2+1))

    Deformation = set(Deformation)
    DeformationFinal = []
    for p in Deformation:
        if p in PermutationDict:
            # DeformationFinal+=list(PermutationDict[p])
            del PermutationDict[p]
            DeformationFinal.append(p)

    print "remaining length of permutation dictionary:", len(PermutationDict)
    return list(DeformationFinal)

def get_Unique_Permutation(permutationList, TimeRotation=True):
    Order = len(permutationList[0])/2
    PermutationDict={}
    for p in permutationList:
        PermutationDict[tuple(p)]=None
    for per in permutationList:
        if not PermutationDict.has_key(tuple(per)):
            continue
        Deformation = [per]

        if TimeRotation:
            for idx in range(1, Order):
                for i in range(len(Deformation)):
                    for j in range(1, idx):
                        Deformation.append(swap_interaction(Deformation[i], idx*2, idx*2+1, j*2, j*2+1))

        for idx in range(1,Order):
            for i in range(len(Deformation)):
                Deformation.append(swap_LR(Deformation[i], idx*2, idx*2+1))

        # for idx in range(1,Order):
            # for i in range(len(Deformation)):
                # Deformation.append(swap_LR_Hugen(Deformation[i], idx*2, idx*2+1))

        Deformation = set(Deformation)
        for p in Deformation:
            if tuple(p)==tuple(per):
                continue
            if p in permutationList:
                del PermutationDict[p]

    print "remaining length of permutation dictionary:", len(PermutationDict)
    return PermutationDict.keys()

def Group(PermutationDict, TimeRotation=True):
    UnlabeledDiagramList=[]
    FactorList=[]
    # for permutation in PermutationList[0:1]:
    while len(PermutationDict)>0:
        print "Remaining diagram {0}".format(len(PermutationDict))
        permutation=PermutationDict.keys()[0]
        Deformation=check_Unique_Permutation(permutation, PermutationDict, TimeRotation)
        # if len(Deformation)>0:
        UnlabeledDiagramList.append(Deformation)
    return UnlabeledDiagramList

def FindIndependentK(permutation, reference, InteractionPairs):
    # kList=[(random.randint(0, Nmax), random.randint(0,Nmax)) for i in range(len(InteractionPairs)+1)]
    N=len(InteractionPairs)
    Matrix=np.zeros((2*N,3*N))
    for i in range(2*N):
        interaction=int(i/2)+2*N
        sign=i%2
        Matrix[i,interaction]=-(-1)**sign
        Matrix[i, i]=-1
        Matrix[i, permutation.index(i)]=1
    # print Matrix
    vectors = nullspace(Matrix)
    # print len(vectors)
    # print vectors
    freedoms=vectors.shape[1]
    if freedoms!=N+1:
        print "Warning! Rank is wrong for {0} with \n{1}".format(permutation, vectors)
    return vectors

def AssignMomentums(permutation, reference, InteractionPairs):
    N=len(InteractionPairs)
    vectors=FindIndependentK(permutation, reference, InteractionPairs)
    freedoms=vectors.shape[1]
    karray=np.array([random.random() for _ in range(freedoms)])
    kVector=np.dot(vectors, karray)
    # kVector=vectors[:,0]
    return kVector[:2*N], kVector[2*N:]

def GenerateAllDiagram(UnlabeledDiagram, InteractionPairs, DoesCheck=True):
    AllDiagram=[]
    for d in UnlabeledDiagram:
        Deformation=[d]
        for idx in range(1,Order):
            for i in range(len(Deformation)):
                Deformation.append(swap_LR_Hugen(Deformation[i], idx*2, idx*2+1))

        DeformationFinal=list(Deformation)
        if DoesCheck is True:
            for p in Deformation:
                kG, kW=AssignMomentums(p, Reference, InteractionPairs)

                Flag=True
                for i in range(len(kW)):
                    if Flag and abs(kW[i])<1e-12:
                        # print "k=0 on W {0}: {1}".format(p, kW[i])
                        DeformationFinal.remove(p)
                        Flag=False
                        break

                for j in range(1,len(kW)):
                    if Flag and abs(abs(kW[0])-abs(kW[j]))<1e-12:
                        # start=2*i
                        # end=p[p[start]]
                        # if start==end and i!=0:
                            # continue
                        # start=2*i+1
                        # end=p[p[start]]
                        # if start==end and i!=0:
                            # continue
                        # print "Same k on W for {0}: {1} on {2}; {3} on {4}".format(p, kW[i],i,kW[j],j)
                        DeformationFinal.remove(p)
                        Flag=False
                        break

                # for i in range(0,len(kG)):
                    # for j in range(i+1,len(kG)):
                        # if Flag and abs(kG[i]-kG[j])<1e-12:
                            # # print "Same k on G for {0}: {1} on {2}; {3} on {4}".format(p, kG[i],i,kG[j],j)
                            # # print "Same k on W for {0}: {1}; 1, {2}".format(p, kG[i],kG[j])
                            # DeformationFinal.remove(p)
                            # Flag=False
                            # # print "Flag",Flag
                            # break
        AllDiagram+=DeformationFinal
    return AllDiagram

def GenerateAllFreeEnergyDiagram(UnlabeledDiagram, InteractionPairs):
    # print "Order", Order
    # print "Diagram", UnlabeledDiagram
    d=UnlabeledDiagram
    Deformation=[d]
    for idx in range(Order):
        for i in range(len(Deformation)):
            Deformation.append(swap_LR_Hugen(Deformation[i], idx*2, idx*2+1))

    for idx in range(Order):
        for i in range(len(Deformation)):
            for j in range(idx):
                Deformation.append(swap_interaction(Deformation[i], idx*2, idx*2+1, j*2, j*2+1))

    for idx in range(Order):
        for i in range(len(Deformation)):
            Deformation.append(swap_LR(Deformation[i], idx*2, idx*2+1))
    DeformationFinal=[tuple(e) for e in Deformation]
    return DeformationFinal

def StartPoint(Order):
    StartPoint=range(Order*2)
    Momentum=np.zeros([Order+1, 2*Order], dtype=int)
    Momentum[0,0]=1
    Momentum[-1,-1]=1
    for i in range(1,Order):
        StartPoint[i*2-1], StartPoint[i*2]=StartPoint[i*2], StartPoint[i*2-1]
        Momentum[i, i*2-1]=1
        Momentum[i, i*2]=1
    FermiSign=-1  #n+1 loop  contributes (-1)^(n+1) and order n contributes (-1)^n 
    return tuple(StartPoint), Momentum, FermiSign

def CheckConservation(permutation, MomentumBases, IsPolarization=False):
    Order=len(permutation)/2
    if matrix_rank(MomentumBases)!=Order+1:
        print "rank is wrong with permutation {0}\n{1}".format(permutation, MomentumBases)
        return False
    Momentum=np.zeros(2*Order)
    for i in range(Order):
        Momentum+=random.random()*MomentumBases[i,:]
    # print len(Momentum)
    for i in range(Order):
        In1, In2=2*i, 2*i+1
        Out1=permutation.index(2*i)
        Out2=permutation.index(2*i+1)
        TotalMom=Momentum[In1]+Momentum[In2]-Momentum[Out1]-Momentum[Out2]
        if abs(TotalMom)>1e-10:
            print "Vertex {0} breaks the conservation laws. Bases: \n{1}".format(i, Momentum)
            print In1, In2, Out1, Out2
            print permutation
            print MomentumBases
            return False
    if IsPolarization:
        #the first loop basis has to be the external momentum
        Ext=np.zeros(Order+1, dtype=int)
        Ext[0]=1
        if not np.all(MomentumBases[:,0]-MomentumBases[:,permutation.index(0)]-Ext==0):
            print "The first loop basis is not the external momentum"
            print permutation, MomentumBases
            sys.exit(0)
        if not np.all(MomentumBases[:,1]-MomentumBases[:,permutation.index(1)]+Ext==0):
            print "The first loop basis is not the external momentum"
            print permutation, MomentumBases
            sys.exit(0)
    return True

def GenerateMomentum(permutation, OldMomentum, i, j):
    if i/2==j/2:
        return None
    Order=len(permutation)/2
    if permutation[i]/2==permutation[j]/2:
        Momentum=np.copy(OldMomentum)
    else:
        Momentum=np.copy(OldMomentum)
        ni=i/2
        nj=j/2
        ip=4*ni+1-i
        jp=4*nj+1-j
        Momentum[:,[i,j]]=Momentum[:,[j,i]]

        if permutation[ip]==jp or permutation[ip]==j:
            Momentum[:,ip]+=Momentum[:,j]-Momentum[:,i]
            # print "Connect ip to jp", i, j, ip, jp, permutation[ip], permutation[jp]
        elif permutation[jp]==i or permutation[jp]==ip:
            Momentum[:,jp]+=Momentum[:,i]-Momentum[:,j]
            # print "Connect jp to ip", i, j, ip, jp, permutation[ip], permutation[jp]
        else:
            return -1
    if not CheckConservation(permutation, Momentum):
        print "Conservation or Rank Check fails."
        sys.exit(0)
        return None
    return Momentum

def GetAllPermutations(Order, DiagDict, DiagInvDict, DiagSymDict):
    """ 
    output:
        Diagrams: a dictionary contains a map from the original diagram to the diagram with optimized bases
        OptDiagrams: a dictionary contains a map from the optimized diagram to a list (original diagram, momentum bases for the optimized diagram, the symmetry factor for the optimized diagram)
    """
    reference=GetReference(Order)
    InteractionPairs=GetInteractionPairs(Order)

    permutation, Momentum, FermiSign=StartPoint(Order)
    PermuList = [permutation]
    MomList=[Momentum]
    SignList=[FermiSign]

    idx = 0
    while idx < 2*Order:
        print "Index {0}".format(idx)
        for i in range(len(PermuList)):
            for j in range(idx):
                newpermutation=tuple(swap(PermuList[i], idx, j))
                # print "Old {0}, New {1} by switching {3},{4}\n OldBases: \n{2}".format(PermuList[i], newpermutation, MomList[i], idx, j)
                if IsConnected(newpermutation, reference, InteractionPairs):
                    Momentum=GenerateMomentum(newpermutation, MomList[i], idx, j) 
                    if Momentum is not None and Momentum is not -1:
                        # print "old :{0}, old_bases: \n{1}\n new:{2}, switch {4},{5},  new_bases: \n {3}".format(PermuList[i], MomList[i], newpermutation, Momentum, idx, j)
                        PermuList.append(newpermutation)
                        MomList.append(Momentum)
                        SignList.append(-SignList[i])
                        # print newpermutation, -SignList[i]
                    # else:
                        # if FreeEnergyDiagramInvDict[newpermutation]==(2,5,4,7,6,9,8,1,0,3):
                        # PermuList.append(newpermutation)
                        # MomList.append(None)
                        # SignList.append(-SignList[i])
                        # print "Can not generate Momentum"
                            # print "old: ", PermuList[i], "new", newpermutation
                            # print "Permutate", idx, j
                            # print "OldBases:\n", MomList[i]
                        # sys.exit(0)
        idx += 1
    
    Diagrams={}
    for i in range(len(PermuList)):
        p=PermuList[i]
        Diag=tuple(DiagInvDict[p])
        SymFactor=SignList[i]*abs(DiagSymDict[Diag])
        Diagrams[Diag]=(p, MomList[i], SymFactor)

    OptDiagrams={}
    for k in Diagrams.keys():
        print "Diagram {0}: {1} with SymFactor {3}\n {2}".format(k, Diagrams[k][0], Diagrams[k][1], Diagrams[k][2])
        p, Mom, Sym=Diagrams[k]
        OptDiagrams[p]=(k, Mom, Sym)

    print "Total Diagram {0} vs {1}".format(len(Diagrams.keys()),len(DiagDict.keys()))

    for k in DiagDict.keys():
        if Diagrams.has_key(k) is False:
            print k

    return Diagrams, OptDiagrams

def FindAllLoops(permutation):
    order=len(permutation)/2
    Visited=set()
    path=[]
    for e in permutation:
        newloop=[]
        vertex=e
        while vertex not in Visited:
            newloop.append(vertex)
            Visited.add(vertex)
            vertex=permutation[vertex]
        if len(newloop)>0:
            path.append(newloop)
    if sum([len(l) for l in path])!=2*order:
        print "length of all loops should be 2*order"
        sys.exit(0)
    return path

def SymmetrizeLoops(OptDiagrams, DiagDict, DiagInvDict):
    order=len(OptDiagrams.keys()[0])/2
    NewOptDiagrams=dict(OptDiagrams)
    for diag in OptDiagrams.keys():
        if diag not in NewOptDiagrams.keys():
            continue
        Momentum=OptDiagrams[diag][1]
        AllDiag=[list(diag)]
        AllMom=[Momentum]
        Loops=FindAllLoops(diag)
        # print "original",diag
        # print "loops",Loops
        for l in Loops:
            for k in range(len(AllDiag)):
                rDiag=list(AllDiag[k])
                rMom=np.array(AllMom[k])
                for e in l:
                    Next=diag[e]
                    rDiag[Next]=e
                    rMom[:,Next]=-Momentum[:,e]
                # print "generate", rDiag
                AllDiag.append(rDiag)
                AllMom.append(rMom)
        EqualDict={}
        OptDict={}
        for d,m in zip(AllDiag, AllMom):
            if not CheckConservation(d, m):
                print "Conservation check fails!"
                sys.exit(0)
            if DiagInvDict.has_key(tuple(d)):
                if not EqualDict.has_key(DiagInvDict[tuple(d)]):
                    EqualDict[DiagInvDict[tuple(d)]]=(d, m)
        print "{0} equals to {1}".format(diag, [EqualDict[k][0] for k in EqualDict.keys()])

        for k in NewOptDiagrams.keys():
            if k is diag:
                continue
            if EqualDict.has_key(DiagInvDict[tuple(k)]):
                d, m=EqualDict[DiagInvDict[tuple(k)]]
                value=NewOptDiagrams[k]
                new_value=(value[0], m, value[-1])
                NewOptDiagrams[tuple(d)]=new_value
                del NewOptDiagrams[k]

    print "Old", OptDiagrams.keys()
    print "New", NewOptDiagrams.keys()
    return NewOptDiagrams

# def GetWLoopBases(diag, mom):

def SaveToFile(UniqueDiagrams, Name):
    if len(UniqueDiagrams)==0: 
        return
    diag, sym, mom, all_diag=UniqueDiagrams[0]
    order=len(diag)/2
    with open("./Diag{0}{1}.txt".format(Name, order), "a") as f:
        f.write("# DiagNum\n{0}\n\n".format(len(UniqueDiagrams)))
        for diag, sym, mom, all_diag in UniqueDiagrams:
            f.write("# Topology\n")
            for i in diag:
                f.write("{0:2d} ".format(i))
            f.write("\n")
            f.write("# SymmetryFactor\n{0}\n".format(sym))
            f.write("# Loop Bases\n")
            for i in range(order+1):
                for j in range(2*order):
                    f.write("{0:2d} ".format(mom[i,j]))
                f.write("\n")

            print "all diagram to print for ", diag
            # print all_diag

            ###### Generate all Mom ################
            # WMom=[]
            # LoopNum=[]
            # for d in all_diag:
                # # print "iterate", d
                # # print "Basis\n", mom
                # for j in range(1,order):
                    # end=2*j
                    # start=d.index(end)
                    # WMom.append(mom[:,start]-mom[:, end])
                    # # print start, end
                    # # print "mom", mom[:, start]-mom[:, end]

            # for i in range(order+1):
                # for j in range((order-1)*len(all_diag)):
                    # # end=2*j
                    # # start=diag.index(end)
                    # # print start, end, mom[i, start]-mom[i, end]
                    # f.write("{0:2d} ".format(WMom[j][i]))
                # f.write("\n")

            ###### Generate indepdent Mom ################
            WMom=[]
            LoopNum=[]
            AllDiagList=[]
            AllDiagList.append(tuple(diag))
            f.write("#Ver4 Legs: InLeft OutLeft InRight OutRight\n")
            for j in range(1,order):
                end=2*j
                start=diag.index(end)
                WMom.append(mom[:,start]-mom[:, end])
                WMom.append(mom[:,start]-mom[:, end+1])
                TempList=[]
                start1=diag.index(end+1)
                f.write("{0} {1} {2} {3} ".format(start, end, start1, end+1))
                for d in AllDiagList:
                    dtemp=list(d)
                    TempList.append(tuple(d))
                    dtemp[start], dtemp[start1]=dtemp[start1], dtemp[start]
                    TempList.append(tuple(dtemp))
                AllDiagList=TempList
            f.write("\n")

            # print diag
            # print mom
            # print "All diag", AllDiagList
            # print np.array(WMom).T

            f.write("#Ver Loop Bases\n")

            for i in range(order+1):
                for j in range(2*(order-1)):
                    # end=2*j
                    # start=diag.index(end)
                    # print start, end, mom[i, start]-mom[i, end]
                    f.write("{0:2d} ".format(WMom[j][i]))
                f.write("\n")

            f.write("#SpinFactor\n")
            for d in AllDiagList:
                path=FindAllLoops(d)
                nloop=len(path)

                ########### for spin susceptibility   #####################
                print "path", path
                Flag=False
                for p in path:
                    if 0 in p and 1 in p:
                        Flag=True

                if Flag==False:
                    print "false", d, path
                    f.write("{0:2d} ".format(0))
                else:
                    f.write("{0:2d} ".format(-(-2)**nloop))
                #####################################################

                # f.write("{0:2d} ".format(-(-2)**nloop))
                # f.write("{0:2d} ".format(-(-1)**nloop))

            f.write("\n")
            f.write("\n")

        
if __name__=="__main__":
    Order=3

    Order-=1
    Reference=GetReference(Order)
    InteractionPairs=GetInteractionPairs(Order)

    FreeEnergyDiagramDict={}
    FreeEnergyDiagramSymDict={}
    FreeEnergyDiagramInvDict={}
    with open("./Diagram/HugenDiag{0}.txt".format(Order)) as f:
        d=f.read()
        exec(d)
        TotalSym=0
        DiagList=ReNameDiag(Diag, Initial=0)

        print "lnZ diagram List:", DiagList
        print "lnZ diagram Symmetry factor:", Sym

        for d, s in zip(DiagList, Sym):
            AllDiagrams=GenerateAllFreeEnergyDiagram(d, InteractionPairs)
            print "{0} with SymFactor: {1}, and Number: {2} with duplicate {3}".format(d, s, len(AllDiagrams), len(AllDiagrams)/len(set(AllDiagrams))) 
            # print "{0}".format(set(AllDiagrams)) 
            TotalSym+=float(len(AllDiagrams))/abs(s)
            FreeEnergyDiagramDict[tuple(d)]=AllDiagrams
            FreeEnergyDiagramSymDict[tuple(d)]=s
            for e in AllDiagrams:
                FreeEnergyDiagramInvDict[tuple(e)]=tuple(d)
        print "Total Free energy diagrams: {0}, TotalSym: {1}".format(len(FreeEnergyDiagramDict), TotalSym)

        DiagDict, OptDiagDict=GetAllPermutations(Order, FreeEnergyDiagramDict, FreeEnergyDiagramInvDict, FreeEnergyDiagramSymDict)
        # SymmetrizeLoops(OptDiagDict, FreeEnergyDiagramDict, FreeEnergyDiagramInvDict)


        # Order=6
        # permutation=(2,5,4,7,6,9,8,1,0,3)
        # print FreeEnergyDiagramInvDict[permutation]
        # newpermutation=swap(permutation, 0, 6)
        # print "new", newpermutation
        # original=FreeEnergyDiagramInvDict[newpermutation]
        # print "Relevant permuation", DiagDict[original][0]
        # print "MomBases:\n", DiagDict[original][1]

    Order+=1 #get the order for polarization
    Reference=GetReference(Order)
    InteractionPairs=GetInteractionPairs(Order)
    TotalDiagNum=0
    TotalSym=0.0

    TempDiagList=OptDiagDict.keys()
    DiagList=[]
    SymList=[]
    MomList=[]
    for d in TempDiagList:
        DiagList.append(tuple([e+2 for e in d]))
        MomList.append(OptDiagDict[d][1])
        SymList.append(OptDiagDict[d][2])
    print DiagList
    print SymList

    OptPolarDict={}
    PolarDiagramDict={}
    PolarDiagramInvDict={}
    KristjanDiag=[]
    PolarUniqueDiagrams=[]
    os.system("rm DiagPolar{0}.txt".format(Order))

    for d, s, m in zip(DiagList, SymList, MomList): 
        print "working on {0}".format(d)
        # PolarDict=AttachExternalVertex(d)
        PolarDict=AttachExternalVertex(d, m, s)
        # print "PolarDict", PolarDict.keys()

        # print PolarDict

        print "Check Tadpole..."
        for p in PolarDict.keys():
            if HasTadpole(p, Reference):
                del PolarDict[p]

        print "Check Fock..."
        for p in PolarDict.keys():
            if HasFock(p, Reference):
                del PolarDict[p]

        print PolarDict.keys(), len(PolarDict.keys())
        UnlabeledDiagramList =  Group(dict(PolarDict), TimeRotation=False)
        LessPolarDict={}
        for g in UnlabeledDiagramList:
            LessPolarDict[g[0]]=PolarDict[g[0]]
            LessPolarDict[g[0]][1]/=len(g)
        UnlabeledDiagramList =  Group(dict(LessPolarDict), TimeRotation=True)

        print "Total Unique Diagrams for Polar:  {0}\n".format(len(UnlabeledDiagramList))
        print UnlabeledDiagramList
        # print "Total Unique Buble Diagrams for Polar: {0}\n".format(len(UnlabeledBubleDiagramList))

        # NewUnlabeledDiagram=[]
        # for g in UnlabeledDiagramList:
            # EqualDict={}
            # for e in g:
                # TimeSeries=tuple(LessPolarDict[e][-1])
                # if not EqualDict.has_key(TimeSeries):
                    # EqualDict[TimeSeries]=[e, 1]
                # else:
                    # EqualDict[TimeSeries][-1]+=1
                    # del LessPolarDict[e]
            # diagrams=[]
            # for k in EqualDict.keys():
                # d=EqualDict[k][0]
                # diagrams.append(d)
                # LessPolarDict[d][1]/=EqualDict[k][1]
            # NewUnlabeledDiagram.append(diagrams)
        # UnlabeledDiagramList=NewUnlabeledDiagram

        TempDiagList=[]
        NumList=[]
        indexList=[]
        index=0
        for g in UnlabeledDiagramList:
            indexList.append(range(index, index+len(g)))
            index+=len(g)
            for e in g:
                TempDiagList.append(e)
                NumList.append(len(g))

        coeffMatrix=np.zeros([len(TempDiagList), len(TempDiagList)])
        coeffMatrixN=np.zeros([len(TempDiagList), len(TempDiagList)])
        coeffdict={}
        for i in range(len(TempDiagList)):
            for j in range(len(TempDiagList)):
                # if i==j:
                    # continue
                di=tuple(TempDiagList[i])
                dj=tuple(TempDiagList[j])
                si=LessPolarDict[di][1]
                sj=LessPolarDict[dj][1]
                val=sum(LessPolarDict[di][-1]*LessPolarDict[dj][-1]/float(si)/float(sj))/DUPLICATE
                coeffMatrix[i][j]=val
                coeffMatrixN[i][j]=val*NumList[i]*NumList[j]
        # print "Coeff Matrix:\n", coeffMatrixN
        # print "Sum: ", np.sum(coeffMatrix)
        # Maximum=10000.0
        # Best=None
        # for m in itertools.product(*indexList):
            # val=np.sum(np.sum(coeffMatrixN[np.ix_(m,m)]))
            # if val<Maximum:
                # Best=m
                # Maximum=val
        # print "Best: ", m, Maximum

        UniqueDiagrams=[]
        IrreducibleDiagrams=[]

        index=0
        for g in UnlabeledDiagramList:
            for e in g:
                print "{0}: {0}".format(e, coeffMatrixN[index, :])
                index+=1
                # UniqueDiagrams.append(e)
            print "Total {0} with symmetry factor {1}\n".format(len(g), float(1.0)/LessPolarDict[g[0]][1])
            # TotalSym+=abs(float(len(g))/s)
            #count how many diagrams
            AllDiagrams=GenerateAllDiagram([g[0], ], InteractionPairs, True)

            print "all diagram {0}:{1}".format(g[0], AllDiagrams)
            if len(AllDiagrams)==0:
                print "Self-energy", g[0]
                continue
            TotalDiagNum+=1
            UniqueDiagrams.append((g[0], float(len(g))/LessPolarDict[g[0]][1], LessPolarDict[g[0]][0], AllDiagrams))
            PolarUniqueDiagrams.append((g[0], float(len(g))/LessPolarDict[g[0]][1], LessPolarDict[g[0]][0], AllDiagrams))

            KristjanDiag+=get_Unique_Permutation(AllDiagrams)

            #generate all polarization diagrams
            # TotalSym+=len(AllDiagrams)*abs(float(len(g))/s)
            TotalSym+=len(AllDiagrams)*abs(float(len(g))/LessPolarDict[g[0]][1])
            AllDiagrams=GenerateAllDiagram(AllDiagrams, InteractionPairs, False)


            PolarDiagramDict[tuple(g[0])]=AllDiagrams
            for d in AllDiagrams:
                PolarDiagramInvDict[tuple(d)]=tuple(g[0])
            OptPolarDict[tuple(g[0])]=(d,PolarDict[tuple(g[0])][0],s)

        # SymmetrizeLoops(OptPolarDict, PolarDiagramDict, PolarDiagramInvDict)
            # print "All diagrams: ", AllDiagrams


    SaveToFile(PolarUniqueDiagrams, "Polar")


    print "Total: {0}, TotalSym: {1}".format(TotalDiagNum, TotalSym)
    index=0
    with open("./Diag{0}.txt".format(Order), "w") as f:
        KristjanDiag.sort()
        print len(KristjanDiag)
        for g in KristjanDiag:
            # print index, g
            f.write("{0} {1}\n".format(index, g))
            index+=1



    # sys.exit(0)







        

        








