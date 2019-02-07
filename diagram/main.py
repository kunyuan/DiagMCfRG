from free_energy import *
from polar import *

if __name__ == "__main__":

    Order = 4

    LnZOrder = Order-1
    DiagFile = "./Diagram/HugenDiag{0}.txt".format(LnZOrder)
    LnZ = free_energy(LnZOrder)
    # Load pre-generated lnZ diagrams
    # build labeled Feynman diagram to unlabled Hugenholtz diagram mapping
    LnZ.LoadDiagrams(DiagFile)

    OptLnZHugenDiagList = LnZ.OptimizeLoopBasis()

    Polar = polar(Order)

    OptPolarHugenDiagList = [Polar.AttachExtVer(
        d) for d in OptLnZHugenDiagList]
