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
    int DiagNum;
    return Ver4Loop(InL, InR, DirTran, DirTran, LoopNum, 0, 3,
                    1,  // calculate u diagram only
                    -1, // do RG diagram
                    Level, DiagNum);
  }

  //   momentum Internal = Var.LoopMom[3];
  //   momentum Internal2 = Var.LoopMom[3] + DirTran;
  //   const momentum &VerLExTran = InL - Internal - DirTran;

  //   const momentum &VerRInL = Internal2;
  //   const momentum &VerRExTran = Internal - InR;

  //   double Weight = 0.0;
  //   double VerWeight = 0.0;
  //   double GWeight = 0.0;
  //   double TauR2L, TauL2R;
  //   int DiagIndex = 0;
  //   int DiagNum = 0;
  //   int Level = 0;
  //   for (int loop = 0; loop < LoopNum; loop++) {
  //     int LTauIndex = 0;
  //     int RTauIndex = LoopNum + 1 - (loop + 1);

  //     //====================  DIRECT  Diagram =============================
  //     // left vertex
  //     int LIndex = DiagNum;
  //     Ver4Loop(InL, Internal, DirTran, VerLExTran, loop, LTauIndex, 4, Level,
  //              DiagNum);
  //     int LDiagNum = DiagNum - LIndex;

  //     // right vertex
  //     int RIndex = DiagNum;
  //     Ver4Loop(VerRInL, InR, DirTran, VerRExTran, LoopNum - 1 - loop,
  //     RTauIndex,
  //              4 + loop, Level, DiagNum);
  //     int RDiagNum = DiagNum - RIndex;

  //     for (int left = LIndex; left < LIndex + LDiagNum; left++)
  //       for (int right = RIndex; right < RIndex + RDiagNum; right++) {
  //         VerWeight =
  //             (_Weight[Level][left][DIRECT] - _Weight[Level][left][EXCHANGE])
  //             *
  //             (_Weight[Level][right][DIRECT] -
  //             _Weight[Level][right][EXCHANGE]);

  //         TauR2L = _ExtTau[Level][left][INR] - _ExtTau[Level][right][OUTL];
  //         TauL2R = _ExtTau[Level][right][INL] - _ExtTau[Level][left][OUTR];

  //         GWeight = Fermi.Green(TauL2R, Internal2, UP, 0, Var.CurrScale) *
  //                   Fermi.Green(TauR2L, Internal, UP, 2, Var.CurrScale);
  //         GWeight += Fermi.Green(TauL2R, Internal2, UP, 2, Var.CurrScale) *
  //                    Fermi.Green(TauR2L, Internal, UP, 0, Var.CurrScale);

  //         Weight += GWeight * VerWeight;
  //         // cout << TauR2L << " vs " << TauL2R << endl;
  //         // cout << "deltaT" << Var.Tau[0] << " vs " << Var.Tau[1] << endl;

  //         // cout << VerLWeight[DIRECT] << " ver " << VerRWeight[DIRECT] <<
  //         endl;
  //         // cout << VerWeight << ", g: " << GWeight << ", total: " << Weight
  //         //      << endl;
  //       }
  //     // double VerWeight = VerLWeight[DIRECT] * VerRWeight[DIRECT];
  //   }
  //   Weight *= pow(-1, LoopNum);
  //   return Weight;
}

double weight::Ver4Loop(const momentum &InL, const momentum &InR,
                        const momentum &DirTran, const momentum &ExTran,
                        int LoopNum, int TauIndex, int LoopIndex, int Channel,
                        int Type, int &Level, int &DiagNum) {
  if (LoopNum == 0) {
    return Ver4Loop0(InL, InR, DirTran, ExTran, TauIndex, LoopIndex, Level,
                     DiagNum);
  } else if (LoopNum == 1)
    return Ver4Loop1(InL, InR, DirTran, ExTran, LoopNum, TauIndex, LoopIndex,
                     Channel, Type, Level, DiagNum);
}

double weight::Ver4Loop0(const momentum &InL, const momentum &InR,
                         const momentum &DirTran, const momentum &ExTran,
                         int TauIndex, int LoopIndex, int &Level,
                         int &DiagNum) {
  //   return VerQTheta.Interaction(InL, InR, DirTran, 0, Var.CurrScale) -
  //   VerQTheta.Interaction(InL, InR, ExTran, 0, Var.CurrScale);

  _Weight[Level][DiagNum][DIRECT] =
      VerQTheta.Interaction(InL, InR, DirTran, 0, Var.CurrScale);
  _Weight[Level][DiagNum][EXCHANGE] =
      VerQTheta.Interaction(InL, InR, ExTran, 0, Var.CurrScale);
  _ExtTau[Level][DiagNum][INL] = Var.Tau[TauIndex];
  _ExtTau[Level][DiagNum][OUTL] = Var.Tau[TauIndex];
  _ExtTau[Level][DiagNum][INR] = Var.Tau[TauIndex];
  _ExtTau[Level][DiagNum][OUTR] = Var.Tau[TauIndex];
  DiagNum += 1;
  return _Weight[Level][DiagNum][DIRECT] - _Weight[Level][DiagNum][EXCHANGE];
}

double weight::Ver4Loop1(const momentum &InL, const momentum &InR,
                         const momentum &DirTran, const momentum &ExTran,
                         int LoopNum, int TauIndex, int LoopIndex, int Channel,
                         int Type, int &Level, int &DiagNum) {

  momentum Internal = Var.LoopMom[LoopIndex];
  double Weight = 0.0;
  double VerWeight = 0.0;
  double LVerWeight, RVerWeight;
  double GWeight = 0.0;
  double TauR2L, TauL2R;
  int DiagIndex = 0;

  // u diagram
  if (Channel == 0 || Channel == 1) {
    momentum Internal2 = Internal + DirTran;
    const momentum &VerLExTran = InL - Internal - DirTran;

    const momentum &VerRInL = Internal2;
    const momentum &VerRExTran = Internal - InR;

    for (int loop = 0; loop < LoopNum; loop++) {
      int LTauIndex = TauIndex;
      int RTauIndex = TauIndex + LoopNum + 1 - (loop + 1);

      //====================  DIRECT  Diagram =============================
      // left vertex
      int LIndex = DiagNum;
      Ver4Loop(InL, Internal, DirTran, VerLExTran, loop, LTauIndex,
               LoopIndex + 1,
               0, // calculate u, s, t sub-ver-diagram
               0, // type 0
               Level, DiagNum);
      int LDiagNum = DiagNum - LIndex;

      // right vertex
      int RIndex = DiagNum;
      Ver4Loop(VerRInL, InR, DirTran, VerRExTran, LoopNum - 1 - loop, RTauIndex,
               LoopIndex + 1 + loop,
               0, // calculate u, s, t sub-ver-diagram
               0, // type 0
               Level, DiagNum);
      int RDiagNum = DiagNum - RIndex;

      for (int left = LIndex; left < LIndex + LDiagNum; left++) {
        for (int right = RIndex; right < RIndex + RDiagNum; right++) {
          VerWeight =
              (_Weight[Level][left][DIRECT] - _Weight[Level][left][EXCHANGE]) *
              (_Weight[Level][right][DIRECT] - _Weight[Level][right][EXCHANGE]);

          TauR2L = _ExtTau[Level][left][INR] - _ExtTau[Level][right][OUTL];
          TauL2R = _ExtTau[Level][right][INL] - _ExtTau[Level][left][OUTR];

          GWeight = Fermi.Green(TauL2R, Internal2, UP, 0, Var.CurrScale) *
                    Fermi.Green(TauR2L, Internal, UP, 2, Var.CurrScale);
          GWeight += Fermi.Green(TauL2R, Internal2, UP, 2, Var.CurrScale) *
                     Fermi.Green(TauR2L, Internal, UP, 0, Var.CurrScale);

          Weight += GWeight * VerWeight * pow(-1, LoopNum);
          // cout << TauR2L << " vs " << TauL2R << endl;
          // cout << "deltaT" << Var.Tau[0] << " vs " << Var.Tau[1] << endl;

          // cout << VerLWeight[DIRECT] << " ver " << VerRWeight[DIRECT] <<
          // endl; cout << VerWeight << ", g: " << GWeight << ", total: " <<
          // Weight
          //      << endl;
        }
        // double VerWeight = VerLWeight[DIRECT] * VerRWeight[DIRECT];
      }
    }
  }
  //   Weight *= pow(-1, LoopNum);
  return Weight;
}