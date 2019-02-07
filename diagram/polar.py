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

    def __GetInteractionPairs(self):
        return tuple([(2*i, 2*i+1) for i in range(1, self.Ver4Num)])

    def __GetReference(self):
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
        ZMomentum = FreeEnergyDiag.Momentum
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

            Assert(diag.CheckConservation(d, Momentum, self.__GetInteractionPairs(
            )), "Momentum does not conserve or rank is wrong!")

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

                if Permutation[1] != Index and Permutation[1] != d.Mirror(Index):
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

            Assert(diag.CheckConservation(Permutation, Mom, self.__GetInteractionPairs(
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

        OptPolarDiagList = []

        for p in PolarDict.keys():
            d = self.BuildADiagram()

            d.Permutation = p
            d.LoopBasis = PolarDict[p][0]
            d.SymFactor = PolarDict[p][1]
            d.VerBasis = [self.__GetReference(), d.GetPermu()]

            OptPolarDiagList.append(d)

    def Group(self, PermutationDict, TimeRotation=True):
        """find the topogically same diagrams in the dictionary"""
        UnlabeledDiagramList = []
        # for permutation in PermutationList[0:1]:
        while len(PermutationDict) > 0:
            print "Remaining diagram {0}".format(len(PermutationDict))
            permutation = PermutationDict.keys()[0]
            Deformation = self.check_Unique_Permutation(
                permutation, PermutationDict, TimeRotation)
            # if len(Deformation)>0:
            UnlabeledDiagramList.append(Deformation)
        return UnlabeledDiagramList

    def check_Unique_Permutation(self, permutation, PermutationDict, TimeRotation):
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

        print "remaining length of permutation dictionary:", len(
            PermutationDict)
        return list(DeformationFinal)

    def get_Unique_Permutation(self, permutationList, TimeRotation=True):
        Order = self.Order
        PermutationDict = {}
        for p in permutationList:
            PermutationDict[tuple(p)] = None
        for per in permutationList:
            if not PermutationDict.has_key(tuple(per)):
                continue
            Deformation = [per]

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

            # for idx in range(1,Order):
                # for i in range(len(Deformation)):
                    # Deformation.append(swap_LR_Hugen(Deformation[i], idx*2, idx*2+1))

            Deformation = set(Deformation)
            for p in Deformation:
                if tuple(p) == tuple(per):
                    continue
                if p in permutationList:
                    del PermutationDict[p]

        print "remaining length of permutation dictionary:", len(
            PermutationDict)
        return PermutationDict.keys()
