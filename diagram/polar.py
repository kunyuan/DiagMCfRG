import diagram as diag
import numpy as np
from logger import *


class polar():
    def __init__(self, Order):
        self.Order = Order
        self.GNum = 2*self.Order
        self.Ver4Num = self.Order
        self.VerNum = 2*self.Ver4Num

        self.ExtLegNum = 2
        self.ExtLoopNum = 1

        self.LoopNum = self.Order+self.ExtLoopNum

    def GetInteractionPairs(self):
        return tuple([(2*i, 2*i+1) for i in range(1, self.Ver4Num)])

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

        Title = "#DiagNum: {0}\n".format(len(PolarHugenList))
        Title += "#Order: {0}\n".format(self.Order)
        Title += "#LoopNum: {0}\n".format(self.LoopNum)
        Title += "#ExtLoopNum: {0}\n".format(self.ExtLoopNum)
        Title += "#Type: {0}\n".format("Normal")
        Title += "\n"

        Body = ""
        for Diag in PolarHugenList:
            Permutation = Diag.GetPermu()
            Mom = Diag.LoopBasis

            Body += "#Topology\n"
            for i in Permutation:
                Body += "{0:2d} ".format(i)
            Body += "\n"

            Body += "#Propagator Type\n"
            for i in Permutation:
                Body += "{0:2d} ".format(0)
            Body += "\n"

            Body += "#Symmetry Factor\n{0}\n".format(Diag.SymFactor)

            Body += "# Loop Basis\n"
            for i in range(self.LoopNum):
                for j in range(self.GNum):
                    Body += "{0:2d} ".format(Diag.LoopBasis[i, j])
                Body += "\n"

            Body += "#Ver4 Legs: InLeft OutLeft InRight OutRight\n"
            for j in range(1, self.Order):
                end1, end2 = 2*j, 2*j+1
                start1 = Permutation.index(end1)
                start2 = Permutation.index(end2)
                Body += "{0} {1} {2} {3} ".format(start1, end1, start2, end2)
            Body += "\n"

            Body += "# Interaction Type\n"
            for i in range(2*(self.Order-1)):
                Body += "{0:2d} ".format(0)
            Body += "\n"

            Body += "#Ver Loop Basis\n"
            InteractionMom = []
            for j in range(1, self.Order):
                end1, end2 = 2*j, 2*j+1
                start1 = Permutation.index(end1)
                start2 = Permutation.index(end2)
                InteractionMom.append(Mom[:, start1]-Mom[:, end1])
                InteractionMom.append(Mom[:, start1]-Mom[:, end2])

            for i in range(self.LoopNum):
                for j in range(2*(self.Order-1)):
                    Body += "{0:2d} ".format(InteractionMom[j][i])
                Body += "\n"

            Body += "#SpinFactor\n"
            for FeynPermu in self.HugenToFeyn(Diag.GetPermu()):
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
                #     Body += "{0:2d} ".format(-(-2)**nloop*(-1)**self.Order)
                #####################################################

                Body += "{0:2d} ".format(2**nloop*int(Sign))
            #   Body += "{0:2d} ".format(-(-1)**nloop)

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
