from free_energy import *

if __name__ == "__main__":

    Order = 3

    DiagFile = "./Diagram/HugenDiag{0}.txt".format(Order)
    LnZ = free_energy(Order)
    # Load pre-generated lnZ diagrams
    # build labeled Feynman diagram to unlabled Hugenholtz diagram mapping
    LnZ.BuildDiagrams(DiagFile)

    LnZ.BuildLoopBasis()
