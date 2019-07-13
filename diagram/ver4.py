import diagram as diag
import numpy as np
from logger import *
import unionfind


class ver4():
    def __init__(self, Order):
        self.Order = Order
        self.GNum = 2*self.Order
        self.Ver4Num = self.Order-1
        self.VerNum = 2*self.Ver4Num

        self.ExtLegNum = 2
        self.ExtLoopNum = 1

        self.LoopNum = self.Order+self.ExtLoopNum

    def GetInteractionPairs(self, WithMeasuring=False):
        if WithMeasuring:
            return tuple([(2*i, 2*i+1) for i in range(0, self.Ver4Num+1)])
        else:
            return tuple([(2*i, 2*i+1) for i in range(1, self.Ver4Num+1)])

    def GetReference(self):
        return tuple(range(self.GNum))

    def BuildADiagram(self):
        d = diag.diagram(self.Order)
        d.Type = "Polar"
        d.GNum = self.GNum
        d.Ver4Num = self.Ver4Num
        d.VerNum = self.VerNum
        d.LoopNum = self.LoopNum
        d.ExtLeg = [0, 1]
        d.ExtLegNum = 2
        d.ExtLoop = [0]
        d.ExtLoopNum = 1
        return d

    def AttachExtVer(self, FreeEnergyDiag):
        PolarDict = dict()
        SHIFT = 2
        Assert(FreeEnergyDiag.Order+1 == self.Order,
               "Polarization and Free energy order should match! {0} vs {1}".format(self.Order, FreeEnergyDiag.Order))

        Diag = FreeEnergyDiag.GetPermu()
        # shift all vertex index with 2
        Diag = [e+SHIFT for e in Diag]

        ZMomentum = FreeEnergyDiag.LoopBasis
        Sym = FreeEnergyDiag.SymFactor

        for i in range(2, len(Diag)+2):
            # Initialization
            # d[i]<== 1 <== 0 <== i
            d = [0, 1]+list(Diag)
            d[1] = d[i]
            d[0] = 1
            d[i] = 0

            Momentum = np.zeros([self.LoopNum, self.GNum], dtype=int)
            Momentum[1:, 2:] = ZMomentum
            Momentum[1:, 0] = ZMomentum[:, i-SHIFT]
            Momentum[1:, 1] = ZMomentum[:, i-SHIFT]
            Momentum[0, 0] = 1

            Assert(diag.CheckConservation(d, Momentum, self.GetInteractionPairs(
            )), "For the first diagram, Momentum does not conserve or rank is wrong!")

            # print "Start with: ", d
            PolarDict[tuple(d)] = [Momentum, Sym]
            ToVisit = [d[1], diag.Mirror(d[1])]
            StartPermu = [tuple(d), tuple(d)]
            StartMom = [Momentum, Momentum]
            Visited = [0]
            # print StartMom[0]
            while len(ToVisit) != 0:
                Index = ToVisit.pop()
                Permutation = list(StartPermu.pop())
                Mom = np.copy(StartMom.pop())
                if Index in Visited:
                    continue

                if Permutation[1] != Index and Permutation[1] != diag.Mirror(Index):
                    print "wrong!", Permutation, Index
                    sys.exit()
                # NextVertex<==1<===PreVertex, Target<======Index
                # NextVertex<=======PreVertex, Target<==1<==Index,
                Target = Permutation[Index]
                NextVertex = Permutation[1]
                PrevVertex = Permutation.index(1)
                Permutation[1] = Target
                Permutation[PrevVertex] = NextVertex
                Permutation[Index] = 1

                deltaMom = np.copy(Mom[:, PrevVertex]-Mom[:, 1])
                Mom[:, 1] = Mom[:, Index]
                Mom[:, Index] += deltaMom

                Assert(diag.CheckConservation(Permutation, Mom, self.GetInteractionPairs(
                )), "Momentum does not conserve or rank is wrong!")

                PolarDict[tuple(Permutation)] = [Mom, Sym]

                Visited.append(Index)

                if Target not in Visited:
                    ToVisit.append(Target)
                    ToVisit.append(diag.Mirror(Target))
                    StartPermu.append(tuple(Permutation))
                    StartPermu.append(tuple(Permutation))
                    StartMom.append(Mom)
                    StartMom.append(Mom)
                # print len(Visited)

        OptPolarDiagDict = {}

        for p in PolarDict.keys():
            d = self.BuildADiagram()

            d.Permutation = p
            d.LoopBasis = PolarDict[p][0]
            d.SymFactor = PolarDict[p][1]
            d.VerBasis = [self.GetReference(), d.GetPermu()]

            d.LoopBasis=self.__FixLoopBasis(d.Permutation, d.LoopBasis)

            OptPolarDiagDict[d.GetPermu()] = d

        # print "Find polar", len(PolarDict.keys())
        return OptPolarDiagDict

    def Group(self, PermutationDict, TimeRotation=True):
        """find the topogically same diagrams in the dictionary"""
        PermutationDict = dict(PermutationDict)
        UnlabelDiagDeformList = []
        # for permutation in PermutationList[0:1]:
        while len(PermutationDict) > 0:
            print "Remaining diagram {0}".format(len(PermutationDict))
            permutation = PermutationDict.keys()[0]
            Deformation = self.__FindDeformation(
                permutation, PermutationDict, TimeRotation)
            # if len(Deformation)>0:
            UnlabelDiagDeformList.append(Deformation)
        return UnlabelDiagDeformList

    def __FindDeformation(self, permutation, PermutationDict, TimeRotation):
        Order = self.Order
        Deformation = [permutation]

        if TimeRotation:
            for idx in range(1, Order):
                for i in range(len(Deformation)):
                    for j in range(1, idx):
                        Deformation.append(diag.SwapTwoInteraction(
                            Deformation[i], idx*2, idx*2+1, j*2, j*2+1))

        for idx in range(1, Order):
            for i in range(len(Deformation)):
                Deformation.append(diag.SwapTwoVertex(
                    Deformation[i], idx*2, idx*2+1))

        for idx in range(1, Order):
            for i in range(len(Deformation)):
                Deformation.append(diag.Direct2Exchange(
                    Deformation[i], idx*2, idx*2+1))

        Deformation = set(Deformation)
        DeformationFinal = []
        for p in Deformation:
            if p in PermutationDict:
                # DeformationFinal+=list(PermutationDict[p])
                del PermutationDict[p]
                DeformationFinal.append(p)

        # print "remaining length of permutation dictionary:", len(
        #     PermutationDict)
        return list(DeformationFinal)

    def ToString(self, PolarHugenList):
        if len(PolarHugenList) == 0:
            return

        IrreDiagList = []
        for Diag in PolarHugenList:
            Permutation = Diag.GetPermu()
            Mom = Diag.LoopBasis

            FeynList = self.HugenToFeyn(Diag.GetPermu())
            FactorList = []

            for FeynPermu in FeynList:
                if self.__IsReducibile(FeynPermu, Diag.LoopBasis):
                    FactorList.append(0)
                else:
                    FactorList.append(1)


            if not np.all(np.array(FactorList) == 0):
                print "Irreducible diagram: ", Permutation
                print "including: "
                for i in range(len(FactorList)):
                    if FactorList[i]!=0:
                        print FeynList[i]
            else:
                # print "Reducible diagram: ", Permutation
                continue

            self.__FindVerSubDiag(Diag.GetPermu())

            IrreDiagList.append([Diag, FeynList, FactorList])
        

        print yellow(
            "Irreducible Ver4 Diag Num: {0}".format(len(IrreDiagList)))

        Title = "#Type: {0}\n".format("RG")
        Title += "#DiagNum: {0}\n".format(len(IrreDiagList))
        Title += "#Order: {0}\n".format(self.Order-2)
        Title += "#GNum: {0}\n".format(self.GNum)
        Title += "#Ver4Num: {0}\n".format(self.Ver4Num)
        Title += "#LoopNum: {0}\n".format(self.LoopNum)
        Title += "#ExtLoopIndex: {0} {1} {2}\n".format(0,1,2)
        Title += "#ExtTransferLoopIndex: {0}\n".format(0)
        Title += "#ExtLegLoopIndex: {0} {1}\n".format(1, 2)
        Title += "#DummyLoopIndex: \n"
        Title += "#TauNum: {0}\n".format(self.Ver4Num+2)
        Title += "#ExtTauIndex: {0} {1}\n".format(0, 1)
        Title += "#DummyTauIndex: \n"
        Title += "\n"

        Body = ""
        for Diag, FeynList, FactorList in IrreDiagList:
            Permutation = Diag.GetPermu()
            Mom = Diag.LoopBasis

            print "Save {0}".format(Permutation)

            Body += "# Permutation\n"
            for i in Permutation:
                Body += "{0:2d} ".format(i)
            Body += "\n"

            Body += "# SymFactor\n{0}\n".format(Diag.SymFactor)

            Body += "# GType\n"
            GType=[0,]*len(Permutation)
            GType[0]=-1
            GType[1]=-1
            GType[Permutation.index(0)]=-1
            GType[Permutation.index(1)]=-1
            for i in Permutation:
                Body += "{0:2d} ".format(GType[i])
            Body += "\n"

            Body += "# VertexBasis\n"
            for i in range(self.GNum):
                Body += "{0:2d} ".format(self.__VerBasis(i))
            Body += "\n"
            for i in range(self.GNum):
                Body += "{0:2d} ".format(self.__VerBasis(Permutation[i]))
            Body += "\n"

            Body += "# LoopBasis\n"
            for i in range(self.LoopNum):
                for j in range(self.GNum):
                    Body += "{0:2d} ".format(Diag.LoopBasis[i, j])
                Body += "\n"

            Body += "# Ver4Legs(InL,OutL,InR,OutR)\n"
            for i in range(1, self.Ver4Num+1):
                # skip the external vertexes 0 and 1
                end1, end2 = 2*i, 2*i+1
                start1 = Permutation.index(end1)
                start2 = Permutation.index(end2)
                Body += "{0:2d} {1:2d} {2:2d} {3:2d} |".format(
                    start1, end1, start2, end2)
            Body += "\n"

            # Get interaction Momemtnum list
            InterMom = self.__GetInteractionMom(Permutation, Mom)

            Body += "# WType(Direct,Exchange)\n"
            for i in range(1, self.Ver4Num+1):
                # type1, type2 = 0, 0
                # if np.all(InterMom[2*i-2] == 0):
                #     type1 = -2
                # if np.all(InterMom[2*i-1] == 0):
                #     type2 = -2
                Body += "{0:2d} {1:2d} |".format(0, 0)
            Body += "\n"

            Body += "# SpinFactor\n"

            for idx, FeynPermu in enumerate(FeynList):
                Path = diag.FindAllLoops(FeynPermu)
                nloop = len(Path)

                Sign = (-1)**nloop*(-1)**(self.Order-1) / \
                    (Diag.SymFactor/abs(Diag.SymFactor))
                # make sure the sign of the Spin factor of the first diagram is positive

                ########### for spin susceptibility   #####################
                # Flag = False
                # for p in Path:
                #     if 0 in p and 1 in p:
                #         Flag = True

                # if Flag == False:
                #     Body += "{0:2d} ".format(0)
                # else:
                #     Body += "{0:2d} ".format(2**nloop *
                #                              int(Sign)*FactorList[idx])
                #####################################################

                # Body += "{0:2d} ".format(2**nloop*int(Sign)*FactorList[idx])
                Body += "{0:2d} ".format(int(Sign)*FactorList[idx])
            #   Body += "{0:2d} ".format(-(-1)**nloop*Factor)

            Body += "\n"
            Body += "\n"

        return Title+Body

    def HugenToFeyn(self, HugenPermu):
        """construct a list of feyn diagram permutation from a hugen diagram permutation"""
        FeynList = []
        FeynList.append(HugenPermu)
        Permutation = HugenPermu
        for j in range(1, self.Order):
            end1, end2 = 2*j, 2*j+1
            start1 = Permutation.index(end1)
            start2 = Permutation.index(end2)

            TempFeynList = []
            for permu in FeynList:
                TempPermu = list(permu)
                TempFeynList.append(tuple(TempPermu))
                TempPermu[start1], TempPermu[start2] = TempPermu[start2], TempPermu[start1]
                TempFeynList.append(tuple(TempPermu))

            FeynList = TempFeynList
        return FeynList

    def __VerBasis(self, index):
        if index <= 1:
            return index
        else:
            return int(index/2)+1
    
    def __IsReducibile(self, Permutation, LoopBasis):
        ExterLoop = [0, ]*self.LoopNum
        ExterLoop[0] = 1
        for i in range(1, self.Ver4Num+1):
            end1, end2 = 2*i, 2*i+1
            start1 = Permutation.index(end1)
            # start2 = Permutation.index(end2)
            VerLoopBasis = LoopBasis[:, start1]-LoopBasis[:, end1]
            # print Permutation, 2*i,  VerLoopBasis

            ####### Check Polarization diagram ##################
            # if(np.array_equal(VerLoopBasis, ExterLoop) or
            #    np.array_equal(-VerLoopBasis, ExterLoop)):
            #     return True

            # remove any hartree insertion
            if(np.all(VerLoopBasis == 0)):
                # print "Contain high-order Hartree: ", Permutation
                return True

        if diag.HasFock(Permutation, self.GetReference()):
            return True

        ###### Check High order Hatree ######################
        kG, kW = diag.AssignMomentums(
            Permutation, self.GetReference(), self.GetInteractionPairs(True))

        # for i in range(len(kW)):
        #     if abs(kW[i]) < 1e-12:
        #             # print "k=0 on W {0}: {1}".format(p, kW[i])
        #         print "Contain high-order Hartree: ", Permutation
        #         return True

        # for j in range(1, len(kW)):
        #     if abs(abs(kW[0])-abs(kW[j])) < 1e-12:
                        # start=2*i
                        # end=p[p[start]]
                        # if start==end and i!=0:
                    # continue
                        # start=2*i+1
                        # end=p[p[start]]
                        # if start==end and i!=0:
                    # continue
                        # print "Same k on W for {0}: {1} on {2}; {3} on {4}".format(p, kW[i],i,kW[j],j)
                # break

        for i in range(0,len(kG)):
            for j in range(i+1,len(kG)):
                if abs(kG[i]-kG[j])<1e-12:
                    # print "Contain Self-energy insertion: ", Permutation
                    # print "Same k on G for {0}: {1} on {2}; {3} on {4}".format(Permutation, kG[i],i,kG[j],j)
                    # print "Same k on W for {0}: {1}; 1, {2}".format(p, kG[i],kG[j])
                    return True
        

        ##### only keep direct diagram #################
        Start=0
        End=Permutation[Start]
        while End!=0 and End!=1:
            Start=End
            End=Permutation[Start]
        if End==1:
            return True

        return False

    def __GetInteractionMom(self, Permutation, Mom):
        InteractionMom = []
        for j in range(1, self.Order):
            end1, end2 = 2*j, 2*j+1
            start1 = Permutation.index(end1)
            start2 = Permutation.index(end2)
            InteractionMom.append(Mom[:, start1]-Mom[:, end1])
            InteractionMom.append(Mom[:, start1]-Mom[:, end2])
        return InteractionMom

        # def get_Unique_Permutation(self, permutationList, TimeRotation=True):
        #     Order = self.Order
        #     PermutationDict = {}
        #     for p in permutationList:
        #         PermutationDict[tuple(p)] = None
        #     for per in permutationList:
        #         if not PermutationDict.has_key(tuple(per)):
        #             continue
        #         Deformation = [per]

        #         if TimeRotation:
        #             for idx in range(1, Order):
        #                 for i in range(len(Deformation)):
        #                     for j in range(1, idx):
        #                         Deformation.append(diag.SwapTwoInteraction(
        #                             Deformation[i], idx*2, idx*2+1, j*2, j*2+1))

        #         for idx in range(1, Order):
        #             for i in range(len(Deformation)):
        #                 Deformation.append(diag.SwapTwoVertex(
        #                     Deformation[i], idx*2, idx*2+1))

        #         # for idx in range(1,Order):
        #             # for i in range(len(Deformation)):
        #                 # Deformation.append(swap_LR_Hugen(Deformation[i], idx*2, idx*2+1))

        #         Deformation = set(Deformation)
        #         for p in Deformation:
        #             if tuple(p) == tuple(per):
        #                 continue
        #             if p in permutationList:
        #                 del PermutationDict[p]

        #     print "remaining length of permutation dictionary:", len(
        #         PermutationDict)
        #     return PermutationDict.keys()

    def __FixLoopBasis(self, Permutation, LoopBasis):
        # if(LoopBasis[1,0]!=1):
        #     if(LoopBasis[1,0]==-1):
        #         LoopBasis[1,:]*=-1
        #     else:
        #         #InL is zero currently
        #         start=2
        #         IsDone=False
        #         while IsDone is False:
                    
        # Path=diag.FindAllLoops(Permutation)
        # length=1000
        # TheOne=None
        # for p in Path:
        #     if(0 in p and len(p)<length):
        #         TheOne=p

        # LoopBasis[0, :]=0
        # for i in TheOne:
        #     LoopBasis[0, i]=1

        return LoopBasis
    
    def __FindVerSubDiag(self, Permutation):
        kG, kW = diag.AssignMomentums(
            Permutation, self.GetReference(), self.GetInteractionPairs(True))
        # InL, OutL, InR, OutR=Permutation[0], Permutation.index(0), Permutation[1], Permutation.index(1)
        SubDiagList=[]
        SubDiagLegList=[]
        SubDiagSizeList=[]

        for i in range(len(Permutation)):
            for j in range(i+1, len(Permutation)):
                for k in range(j+1, len(Permutation)):
                    for l in range(k+1, len(Permutation)):
                        Flag=False
                        if(abs(abs(kG[i]+kG[j])-abs(kG[k]+kG[l])))<1.0e-6:
                            # print "Sub Diagram 1: ", i, j, k, l
                            Flag=True
                        elif(abs(abs(kG[i]-kG[j])-abs(kG[k]+kG[l])))<1.0e-6:
                            # print "Sub Diagram 2: ", i, j, k, l
                            Flag=True
                        elif(abs(abs(kG[i]+kG[j])-abs(kG[k]-kG[l])))<1.0e-6:
                            # print "Sub Diagram 3: ", i, j, k, l
                            Flag=True
                        elif(abs(abs(kG[i]-kG[j])-abs(kG[k]-kG[l])))<1.0e-6:
                            # print "Sub Diagram 4: ", i, j, k, l
                            Flag=True

                        if Flag==True:
                            diagram = set(self.GetInteractionPairs(True))
                            for e in range(len(Permutation)):
                                if(e!=i and e!=j and e!=k and e!=l):
                                    diagram.add((e, Permutation[e]))
                                
                            n_node = len(self.GetInteractionPairs(True))*2

                            diagram_union = unionfind.UnionFind(n_node)

                            for edge in diagram:
                                # print "edge: ", edge
                                if edge[0] != edge[1] and not diagram_union.is_connected(edge[0], edge[1]):
                                    diagram_union.union(edge[0], edge[1])

                            GroupNum=diagram_union.get_n_circles()
                            if GroupNum!=2:
                                Abort("Group number got {0}".format(GroupNum))

                            GroupsMasks=diagram_union.get_circles()
                            ExternIndex=GroupsMasks[0]
                            InternGroup=[node for node in range(len(Permutation)) if GroupsMasks[node]!=ExternIndex]
                            if len(InternGroup)>2:
                                # print "Permu: ", Permutation
                                # print "GroupNum: ", GroupNum
                                # print "Mask:", GroupsMasks
                                Set1=set((Permutation[i], Permutation[j], Permutation[k], Permutation[l]))
                                InL, InR=[Permutation.index(leg) for leg in Set1.intersection(set(InternGroup))]
                                # print "InterGroup:", InternGroup
                                # print "leg", InL, InR
                                # print "\n"
                                End=Permutation[InL]
                                while End in InternGroup:
                                    End=Permutation[End]
                                OutL=Permutation.index(End)

                                End=Permutation[InR]
                                while End in InternGroup:
                                    End=Permutation[End]
                                OutR=Permutation.index(End)

                                Legs=[InL, OutL, InR, OutR]
                                if set(Legs)!=set((i, j, k, l)):
                                    Abort("Legs not equal! {0} vs {1}".format(Legs, (i,j,k,l)))
                                SubDiagLegList.append(Legs)
                                SubDiagList.append(InternGroup)
                                SubDiagSizeList.append(len(InternGroup))

        for index in range(len(SubDiagLegList)):
            print "Subdiagram:", SubDiagList[index]
            print "Legs:", SubDiagLegList[index]
            print "Size:", SubDiagSizeList[index]
                        
    # def __FindDisconnect(self, Permutation, Legs):
    #     start=Legs[0]
    #     Visited=[start,]
    #     ToSearch=[]
    #     if Permutation[start] not in Visited:
    #         ToSearch.append(Permutation[start])
    #     if diag.Mirror(Permutation[start]) not in Visited:
    #         ToSearch.append(diag.Mirror(Permutation[start]))

    #     while len(ToSearch)>0:
    #         start=ToSearch[-1]
    #         Visited.append(start)
    #         del ToSearch[-1]
    #         if start==Legs[1] or start==Legs[2] or start==Legs[3]:
    #             continue
    #         else:
    #             if Permutation[start] not in Visited:
    #                 ToSearch.append(Permutation[start])
    #             if diag.Mirror(Permutation[start]) not in Visited:
    #                 ToSearch.append(diag.Mirror(Permutation[start]))
    #     print "Permutation: ", Permutation
    #     print "Leg: ", Legs
    #     print "Visited: ", Visited




