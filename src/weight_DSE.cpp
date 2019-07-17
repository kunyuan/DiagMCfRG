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
  } else {
    int Level = 0;
    int DiagIndex = 0;
    double Weight = 0.0;
    if (LoopNum == 1) {
      DiagIndex = Ver4Loop(InL, InR, DirTran, LoopNum, 0, 3, DiagIndex, Level,
                           1, // do RG diagram
                           0  // calculate u diagram only
      );
      // Weight = _Weight[Level][0];
    } else if (LoopNum == 2) {
      DiagIndex = Ver4Loop(InL, InR, DirTran, LoopNum, 0, 3, DiagIndex, Level,
                           1, // do RG diagram
                           0, // calculate u diagram only
                           1  // LverOrder
      );
      // Weight = _Weight[Level][2];
      // cout << _Weight[Level][0] << ": " << _Weight[Level][1] << ": "
      //      << _Weight[Level][2] << ": " << _Weight[Level][3] << endl;
    }
    // int count = 0;
    for (int diag = 0; diag < DiagIndex; diag++) {
      Weight += _Weight[Level][diag][1];
      // count++;
    }

    // cout << _Weight[Level][0] << ": " << _Weight[Level][1] << ": "
    //      << _Weight[Level][2] << endl;

    // if (LoopNum == 2)
    //   cout << count << " diag: " << DiagIndex << endl;
    return Weight / pow(40.0, LoopNum) * pow(-1, LoopNum - 1);
  }
}

int weight::Ver4Loop(const momentum &InL, const momentum &InR,
                     const momentum &DirTran, int LoopNum, int TauIndex,
                     int LoopIndex, int DiagIndex, int Level, int Type,
                     int Channel, int LVerOrder) {
  if (LoopNum == 0) {
    return Ver4Loop0(InL, InR, DirTran, TauIndex, LoopIndex, Level, DiagIndex);
  } else if (LoopNum >= 1)
    return Ver4Loop1(InL, InR, DirTran, LoopNum, TauIndex, LoopIndex, DiagIndex,
                     Level, Type, Channel, LVerOrder);
}

int weight::Ver4Loop0(const momentum &InL, const momentum &InR,
                      const momentum &DirTran, int TauIndex, int LoopIndex,
                      int Level, int DiagIndex) {
  //   return VerQTheta.Interaction(InL, InR, DirTran, 0, Var.CurrScale) -
  //   VerQTheta.Interaction(InL, InR, ExTran, 0, Var.CurrScale);
  momentum ExTran = InR + DirTran - InL;
  double DiWeight = VerQTheta.Interaction(InL, InR, DirTran, 0, Var.CurrScale);
  double ExWeight = VerQTheta.Interaction(InL, InR, ExTran, 0, Var.CurrScale);
  _ExtTau[Level][DiagIndex][INL] = Var.Tau[TauIndex];
  _ExtTau[Level][DiagIndex][OUTL] = Var.Tau[TauIndex];
  _ExtTau[Level][DiagIndex][INR] = Var.Tau[TauIndex];
  _ExtTau[Level][DiagIndex][OUTR] = Var.Tau[TauIndex];
  _Weight[Level][DiagIndex][0] = DiWeight - ExWeight;
  _Weight[Level][DiagIndex][1] = 0.0;
  // _Weight[Level][DiagIndex] = DiWeight;

  // ASSERT_ALLWAYS(abs(DirTran.norm() - Var.LoopMom[0].norm()) < 1.0e-5,
  //                "Ext Mom wrong!");

  DiagIndex += 1;
  return DiagIndex;
}

int weight::Ver4Loop1(momentum InL, momentum InR, const momentum &DirTran,
                      int LoopNum, int TauIndex, int LoopIndex, int DiagIndex,
                      int Level, int RG, int Channel, int LVerOrder,
                      int Projection) {

  if (Projection == 1) {
    InL = InL * (Para.Kf / InL.norm());
    InR = InR * (Para.Kf / InR.norm());
    if (RG == 1)
      // derivative of a projection is zero
      return 0;
  }

  momentum Internal = Var.LoopMom[LoopIndex];
  momentum OutL = InL - DirTran;
  momentum OutR = InR + DirTran;
  momentum Internal2, VerLInL, VerLInR, VerLDiTran, VerRInL, VerRInR,
      VerRDiTran;
  double Weight = 0.0;
  double VerWeight = 0.0;
  double LVerWeight, RVerWeight;
  double GWeight = 0.0;
  double GWeightRG = 0.0;
  double TauR2L, TauL2R, SymFactor;
  int nLevel = Level + 1;
  int nDiagIndex = 0;

  for (int chan = 0; chan < 3; chan++) {
    if ((Channel == -1 || Channel == 0) && chan == 0) {
      // t diagram
      Internal2 = Internal + DirTran;
      VerLInL = InL;
      VerLInR = Internal;
      VerLDiTran = DirTran;

      VerRInL = Internal2;
      VerRInR = InR;
      VerRDiTran = DirTran;
      SymFactor = 1.0;
    } else if ((Channel == -1 || Channel == 1) && chan == 1) {
      // u diagram
      Internal2 = Internal - DirTran + InL - InR;
      VerLInL = InL;
      VerLInR = Internal;
      VerLDiTran = InL - OutR;

      VerRInL = Internal2;
      VerRInR = InR;
      VerRDiTran = Internal2 - Internal;
      SymFactor = 1.0;
    } else if ((Channel == -1 || Channel == 2) && chan == 2) {
      // s diagram
      Internal2 = InL + InR - Internal;
      VerLInL = InL;
      VerLInR = InR;
      VerLDiTran = VerLInL - Internal2;

      VerRInL = Internal2;
      VerRInR = Internal;
      VerRDiTran = Internal2 - OutL;
      SymFactor = 0.5;
    } else {
      continue;
    }

    int LVerStart = 0;
    int LVerEnd = LoopNum;
    if (LVerOrder >= 0) {
      LVerStart = LVerOrder;
      LVerEnd = LVerOrder + 1;
    }
    // iterate all possible left vertex order
    for (int loop = LVerStart; loop < LVerEnd; loop++) {
      int LTauIndex = TauIndex;
      int RTauIndex = TauIndex + (loop + 1);

      //====================  DIRECT  Diagram =============================
      // left vertex
      int LIndex = nDiagIndex;
      nDiagIndex = Ver4Loop(VerLInL, VerLInR, VerLDiTran, loop, LTauIndex,
                            LoopIndex + 1, nDiagIndex, nLevel, RG,
                            2 // calculate u, s, t sub-ver-diagram
      );
      int LDiagIndex = nDiagIndex;

      // right vertex
      int RIndex = LDiagIndex;
      nDiagIndex =
          Ver4Loop(VerRInL, VerRInR, VerRDiTran, LoopNum - 1 - loop, RTauIndex,
                   LoopIndex + 1 + loop, nDiagIndex, nLevel, RG,
                   0 // calculate u, s, t sub-ver-diagram
          );
      int RDiagIndex = nDiagIndex;

      for (int left = LIndex; left < LDiagIndex; left++) {
        for (int right = RIndex; right < RDiagIndex; right++) {
          if ((Channel == -1 || Channel == 0) && chan == 0) {
            _ExtTau[Level][DiagIndex][INL] = _ExtTau[nLevel][left][INL];
            _ExtTau[Level][DiagIndex][OUTL] = _ExtTau[nLevel][left][OUTL];
            _ExtTau[Level][DiagIndex][INR] = _ExtTau[nLevel][right][INR];
            _ExtTau[Level][DiagIndex][OUTR] = _ExtTau[nLevel][right][OUTR];
          } else if ((Channel == -1 || Channel == 1) && chan == 1) {
            _ExtTau[Level][DiagIndex][INL] = _ExtTau[nLevel][left][INL];
            _ExtTau[Level][DiagIndex][OUTL] = _ExtTau[nLevel][right][OUTR];
            _ExtTau[Level][DiagIndex][INR] = _ExtTau[nLevel][right][INR];
            _ExtTau[Level][DiagIndex][OUTR] = _ExtTau[nLevel][left][OUTL];
          } else if ((Channel == -1 || Channel == 2) && chan == 2) {
            _ExtTau[Level][DiagIndex][INL] = _ExtTau[nLevel][left][INL];
            _ExtTau[Level][DiagIndex][OUTL] = _ExtTau[nLevel][right][OUTL];
            _ExtTau[Level][DiagIndex][INR] = _ExtTau[nLevel][left][INR];
            _ExtTau[Level][DiagIndex][OUTR] = _ExtTau[nLevel][right][OUTR];
          }
          if (Projection == 1) {
            double avg = (_ExtTau[Level][DiagIndex][INL] +
                          _ExtTau[Level][DiagIndex][OUTL] +
                          _ExtTau[Level][DiagIndex][INR] +
                          _ExtTau[Level][DiagIndex][OUTR]) /
                         4.0;
            _ExtTau[Level][DiagIndex][INL] = avg;
            _ExtTau[Level][DiagIndex][OUTL] = avg;
            _ExtTau[Level][DiagIndex][INR] = avg;
            _ExtTau[Level][DiagIndex][OUTR] = avg;
          }
          if ((Channel == -1 || Channel == 2) && chan == 2) {
            TauR2L = _ExtTau[nLevel][right][INR] - _ExtTau[nLevel][left][OUTR];
            TauL2R = _ExtTau[nLevel][right][INL] - _ExtTau[nLevel][left][OUTL];
          } else {
            TauR2L = _ExtTau[nLevel][left][INR] - _ExtTau[nLevel][right][OUTL];
            TauL2R = _ExtTau[nLevel][right][INL] - _ExtTau[nLevel][left][OUTR];
          }

          _Weight[Level][DiagIndex][0] = 0.0;
          _Weight[Level][DiagIndex][1] = 0.0;

          GWeight = Fermi.Green(TauL2R, Internal2, UP, 0, Var.CurrScale) *
                    Fermi.Green(TauR2L, Internal, UP, 0, Var.CurrScale);
          VerWeight = _Weight[nLevel][left][0] * _Weight[nLevel][right][0];
          _Weight[Level][DiagIndex][0] +=
              GWeight * VerWeight * (-1) * SymFactor;
          Weight += _Weight[Level][DiagIndex][0];

          if (RG == 1) {
            GWeightRG = Fermi.Green(TauL2R, Internal2, UP, 0, Var.CurrScale) *
                        Fermi.Green(TauR2L, Internal, UP, 2, Var.CurrScale);
            GWeightRG += Fermi.Green(TauL2R, Internal2, UP, 2, Var.CurrScale) *
                         Fermi.Green(TauR2L, Internal, UP, 0, Var.CurrScale);

            // if both left and right vertex do not contain derivative, then
            // GG can contain derivative
            VerWeight = _Weight[nLevel][left][0] * _Weight[nLevel][right][0];
            _Weight[Level][DiagIndex][1] +=
                GWeightRG * VerWeight * (-1) * SymFactor;
            Weight += _Weight[Level][DiagIndex][1];

            VerWeight = _Weight[nLevel][left][0] * _Weight[nLevel][right][1];
            VerWeight += _Weight[nLevel][left][1] * _Weight[nLevel][right][0];
            _Weight[Level][DiagIndex][1] +=
                GWeight * VerWeight * (-1) * SymFactor;
            Weight += _Weight[Level][DiagIndex][1];
          }
        }
        DiagIndex++;
      }
    }
  }
  return DiagIndex;
}