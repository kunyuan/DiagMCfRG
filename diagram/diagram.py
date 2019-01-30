import unionfind
from logger import *
from nullspace import rank, nullspace
from numpy.linalg import matrix_rank
import numpy as np


class diagram:
    """a feynman diagram class"""

    def __init__(self, order, permutation=[], legs=[], loop=[], externloop=[], sym=None):
        self.Permutation = permutation
        self.Legs = legs
        self.ExternLoop = externloop
        self.LoopBasis = loop
        self.Order = order
        self.GNum = len(permutation)
        self.Ver4Num = len(permutation)/2
        self.SymFactor = sym

        # reference permutation [0,1,2,3,4,...]
        self.Reference = self.GetReference()
        self.InteractionPairs = self.GetInteractionPairs()

    def IsConnected(self, InteractionPairs):
        diagram = set(self.GetInteractionPairs)
        for i in range(len(self.Permutation)):
            diagram.add((self.Reference[i], self.Permutation[i]))

        n_node = len(InteractionPairs)*2
        diagram_union = unionfind.UnionFind(n_node)

        for edge in diagram:
            if edge[0] != edge[1] and not diagram_union.is_connected(edge[0], edge[1]):
                diagram_union.union(edge[0], edge[1])
        return diagram_union.get_n_circles() == 1

    def GetInteractionPairs(self):
        return [(2*i, 2*i+1) for i in range(self.Ver4Num)]

    def GetReference(self):
        return range(self.GNum)

    def HasTadpole(self):
        for i in range(len(self.Permutation)):
            if self.Reference[i] == self.Permutation[i]:
                return True
        return False

    def HasFock(self):
        for i in range(len(self.Reference)):
            # end=reference[i]
            end = self.Permutation[i]
            if i == 0 or i == 1:
                continue
            if abs(i-end) == 1 and min(i, end) % 2 == 0:
                return True
        return False

    def swap_LR(self, i, j):
        ip, jp = self.Permutation.index(i), self.Permutation.index(j)
        self.Permutation[ip] = j
        self.Permutation[jp] = i
        self.Permutation[i], self.Permutation[j] = self.Permutation[j], self.Permutation[i]

    # def swap_LR_Hugen(permutation, i, j):
    #     permutation = list(permutation)
    #     permutation[i], permutation[j] = permutation[j], permutation[i]
    #     return swap_LR(permutation, i, j)
    #     # return tuple(permutation)

    # def swap_LR_Hugen_Backward(permutation, i, j):
    #     permutation = list(permutation)
    #     permutation[i], permutation[j] = permutation[j], permutation[i]
    #     return permutation
    #     # return tuple(permutation)

    def StartPoint(self):
        self.Permutation = range(self.GNum)
        Momentum = np.zeros([self.GNum/2+1, 2*self.GNum], dtype=int)
        Momentum[0, 0] = 1
        Momentum[-1, -1] = 1
        for i in range(1, self.GNum/2):
            self.Permutation[i*2-1], self.Permutation[i *
                                                      2] = self.Permutation[i*2], self.Permutation[i*2-1]
            Momentum[i, i*2-1] = 1
            Momentum[i, i*2] = 1
        # n+1 loop  contributes (-1)^(n+1) and order n contributes (-1)^n
        FermiSign = -1
        return tuple(StartPoint), Momentum, FermiSign

    def FindAllLoops(self):
        Visited = set()
        path = []
        for e in self.Permutation:
            newloop = []
            vertex = e
            while vertex not in Visited:
                newloop.append(vertex)
                Visited.add(vertex)
                vertex = self.Permutation[vertex]
            if len(newloop) > 0:
                path.append(newloop)
        Assert(sum([len(l) for l in path]) == self.GNum,
               "length of all loops should be 2*order")
        return path
