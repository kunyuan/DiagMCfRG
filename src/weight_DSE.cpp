#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/vector.h"
#include "weight.h"
#include <array>
#include <iostream>
#include <string>

using namespace diag;
using namespace std;

double weight::Evaluate(group &Group) {
  if (Group.Order == 0)
    VerQTheta.Interaction(Var.LoopMom[1], Var.LoopMom[2], Var.LoopMom[0], -1,
                          Var.CurrScale);
  else {
    // return fRG(group & Group);
    // RG loop
  }
}

double weight::fRG(int LoopNum) {

  const momentum &DirTran = Var.LoopMom[0];
  const momentum &InL = Var.LoopMom[1];
  const momentum &InR = Var.LoopMom[2];
  if (LoopNum == 0) {
    // normalization
    return VerQTheta.Interaction(InL, InR, DirTran, -2, Var.CurrScale);
  }

  momentum Internal = Var.LoopMom[3];
  momentum Internal2 = Var.LoopMom[3] + DirTran;
  const momentum &VerLExTran = InL - Internal - DirTran;

  const momentum &VerRInL = Internal2;
  const momentum &VerRExTran = Internal - InR;

  double VerLWeight[2];
  double VerRWeight[2];
  double VerLExtTau[4];
  double VerRExtTau[4];
  double Weight = 0.0;
  for (int loop = 0; loop < LoopNum; loop++) {
    int LTauIndex = 0;
    int RTauIndex = LoopNum + 1 - (loop + 1);

    //====================  DIRECT  Diagram =============================
    Ver4Loop(InL, Internal, DirTran, VerLExTran, loop, LTauIndex, 4, VerLWeight,
             VerLExtTau);
    Ver4Loop(VerRInL, InR, DirTran, VerRExTran, LoopNum - 1 - loop, RTauIndex,
             4 + loop, VerRWeight, VerRExtTau);
    // double
    // VerWeight=(VerLWeight[DIRECT]-VerLWeight[EXCHANGE])*(VerRWeight[DIRECT]-VerRWeight[EXCHANGE])
    double VerWeight = VerLWeight[DIRECT] * VerRWeight[DIRECT];

    double TauR2L = VerLExtTau[INR] - VerRExtTau[OUTL];
    double TauL2R = VerRExtTau[INL] - VerLExtTau[OUTR];

    // cout << TauR2L << " vs " << TauL2R << endl;
    // cout << "deltaT" << Var.Tau[0] << " vs " << Var.Tau[1] << endl;

    double GWeight = Fermi.Green(TauL2R, Internal2, UP, 0, Var.CurrScale) *
                     Fermi.Green(TauR2L, Internal, UP, 2, Var.CurrScale);
    GWeight += Fermi.Green(TauL2R, Internal2, UP, 2, Var.CurrScale) *
               Fermi.Green(TauR2L, Internal, UP, 0, Var.CurrScale);

    Weight += GWeight * VerWeight;
    // cout << VerLWeight[DIRECT] << " ver " << VerRWeight[DIRECT] << endl;
    // cout << VerWeight << ", g: " << GWeight << ", total: " << Weight << endl;
  }
  Weight *= pow(-1, LoopNum);
  return Weight;
}

void weight::Ver4Loop(const momentum &InL, const momentum &InR,
                      const momentum &DirTran, const momentum &ExTran,
                      int LoopNum, int TauIndex, int LoopIndex, double *Weight,
                      double *ExtTau) {
  if (LoopNum == 0) {
    Weight[DIRECT] = VerQTheta.Interaction(InL, InR, DirTran, 0, Var.CurrScale);
    Weight[EXCHANGE] =
        VerQTheta.Interaction(InL, InR, ExTran, 0, Var.CurrScale);
    ExtTau[0] = Var.Tau[TauIndex];
    ExtTau[1] = Var.Tau[TauIndex];
    ExtTau[2] = Var.Tau[TauIndex];
    ExtTau[3] = Var.Tau[TauIndex];
    return;
  }
}