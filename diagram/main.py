from free_energy import *
# from polar import *
from ver4 import *
import copy
import sys

if __name__ == "__main__":
    InputStr = raw_input("Input Diagram Order: ")
    Order = int(InputStr)

    LnZOrder = Order-1
    DiagFile = "./Diagram/HugenDiag{0}.txt".format(LnZOrder)
    LnZ = free_energy(LnZOrder)
    # Load pre-generated lnZ diagrams
    # build labeled Feynman diagram to unlabled Hugenholtz diagram mapping
    print "\nLoad Order {0} LnZ diagrams ...".format(LnZOrder)
    LnZ.LoadDiagrams(DiagFile)

    print red("\nThe optimimal LnZ diagrams:")
    OptLnZHugenDiagList = LnZ.OptimizeLoopBasis()

    Polar = ver4(Order)

    UniqueUnLabelDiagList = []

    for d in OptLnZHugenDiagList:
        print "\n============================================================="
        print blue("Processing LnZ diagram: {0} with SymFactor: {1}".format(
            d.GetPermu(), d.SymFactor))

        print "Attach two external vertexes ..."
        OptPolarHugenDiagDict = Polar.AttachExtVer(d)

        # print "Check Tadpole..."
        # for p in OptPolarHugenDiagDict.keys():
        #     if diag.HasTadpole(p, Polar.GetReference()):
        #         del OptPolarHugenDiagDict[p]

        # print "Check Fock..."
        # for p in OptPolarHugenDiagDict.keys():
        #     if diag.HasFock(p, Polar.GetReference()):
        #         del OptPolarHugenDiagDict[p]

        print "Group polarization diagrams from the same LnZ diagram..."
        UnLabelDiagDeformList = Polar.Group(
            OptPolarHugenDiagDict, TimeRotation=True)
        # each element contains a deforamtion of hugenholz polarization
        # diagrams in the same LnZ group

        print red("Representative polarization Hugenholtz diagram:")
        for d in UnLabelDiagDeformList:
            diagram = copy.deepcopy(OptPolarHugenDiagDict[d[0]])
            diagram.SymFactor = diagram.SymFactor*len(d)
            UniqueUnLabelDiagList.append(diagram)
            print red("{0} with SymFactor {1}".format(
                diagram.GetPermu(), diagram.SymFactor))

    print yellow("Total Unique Polarization diagram: {0}".format(
        len(UniqueUnLabelDiagList)))

    print "Save diagrams ..."
    # with open("./Diag{0}{1}.txt".format("Polar", Order), "w") as f:
    #     f.write(Polar.ToString(UniqueUnLabelDiagList))
    with open("./Diag{0}{1}.txt".format("4Ver", Order-2), "w") as f:
        f.write(Polar.ToString(UniqueUnLabelDiagList))
