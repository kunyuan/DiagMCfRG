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

    print red("\nThe optimimal LnZ diagrams:")
    OptLnZHugenDiagList = LnZ.OptimizeLoopBasis()

    Polar = polar(Order)

    for d in OptLnZHugenDiagList:
        print "\n============================================================="
        print blue("Processing LnZ diagram: {0} with SymFactor: {1}".format(
            d.GetPermu(), d.SymFactor))

        print "Attach two external vertexes ..."
        OptPolarHugenDiagDict = Polar.AttachExtVer(d)

        print "Check Tadpole..."
        for p in OptPolarHugenDiagDict.keys():
            if diag.HasTadpole(p, Polar.GetReference()):
                del OptPolarHugenDiagDict[p]

        print "Check Fock..."
        for p in OptPolarHugenDiagDict.keys():
            if diag.HasFock(p, Polar.GetReference()):
                del OptPolarHugenDiagDict[p]

        # each element contains a deforamtion of hugenholz polarization
        # diagrams in the same LnZ group
        print "Group polarization diagrams from the same LnZ diagram..."
        UnLabeledDiagList = Polar.Group(
            OptPolarHugenDiagDict, TimeRotation=True)
        print red("Representative polarization Hugenholtz diagram:")
        for d in UnLabeledDiagList:
            print red("{0} with SymFactor {1}".format(
                d[0], OptPolarHugenDiagDict[d[0]].SymFactor*len(d)))
