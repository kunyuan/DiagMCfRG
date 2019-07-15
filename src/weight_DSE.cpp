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
    int DiagNum = 0;
    if (LoopNum == 1) {
      Ver4Loop(InL, InR, DirTran, LoopNum, 0, 3,
               0,  // calculate u diagram only
               -1, // do RG diagram
               Level, DiagNum);
    } else if (LoopNum == 2) {
      Ver4Loop(InL, InR, DirTran, LoopNum, 0, 3,
               0,  // calculate u diagram only
               -1, // do RG diagram
               Level, DiagNum, 0);
    }
    double Weight = 0.0;
    for (int diag = 0; diag < DiagNum; diag++)
      if (_Derivative[Level][diag] == true)
        Weight += _Weight[Level][diag];
    return Weight / pow(40.0, LoopNum);
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
                        const momentum &DirTran, int LoopNum, int TauIndex,
                        int LoopIndex, int Channel, int Type, int &Level,
                        int &DiagNum, int LVerOrder) {
  if (LoopNum == 0) {
    return Ver4Loop0(InL, InR, DirTran, TauIndex, LoopIndex, Level, DiagNum);
  } else if (LoopNum >= 1)
    return Ver4Loop1(InL, InR, DirTran, LoopNum, TauIndex, LoopIndex, Channel,
                     Type, Level, DiagNum, LVerOrder);
}

double weight::Ver4Loop0(const momentum &InL, const momentum &InR,
                         const momentum &DirTran, int TauIndex, int LoopIndex,
                         int &Level, int &DiagNum) {
  //   return VerQTheta.Interaction(InL, InR, DirTran, 0, Var.CurrScale) -
  //   VerQTheta.Interaction(InL, InR, ExTran, 0, Var.CurrScale);
  momentum ExTran = InR + DirTran - InL;
  double DiWeight = VerQTheta.Interaction(InL, InR, DirTran, 0, Var.CurrScale);
  double ExWeight = VerQTheta.Interaction(InL, InR, ExTran, 0, Var.CurrScale);
  _ExtTau[Level][DiagNum][INL] = Var.Tau[TauIndex];
  _ExtTau[Level][DiagNum][OUTL] = Var.Tau[TauIndex];
  _ExtTau[Level][DiagNum][INR] = Var.Tau[TauIndex];
  _ExtTau[Level][DiagNum][OUTR] = Var.Tau[TauIndex];
  _Weight[Level][DiagNum] = DiWeight - ExWeight;
  // _Weight[Level][DiagNum] = DiWeight;
  _Derivative[Level][DiagNum] = false;

  // ASSERT_ALLWAYS(abs(DirTran.norm() - Var.LoopMom[0].norm()) < 1.0e-5,
  //                "Ext Mom wrong!");

  DiagNum += 1;
  return _Weight[Level][DiagNum];
}

double weight::Ver4Loop1(const momentum &InL, const momentum &InR,
                         const momentum &DirTran, int LoopNum, int TauIndex,
                         int LoopIndex, int Channel, int Type, int &Level,
                         int &DiagNum, int LVerOrder) {

  momentum Internal = Var.LoopMom[LoopIndex];
  momentum OutL = InL - DirTran;
  momentum OutR = InR + DirTran;
  momentum Internal2, VerLInL, VerLInR, VerLDiTran, VerRInL, VerRInR,
      VerRDiTran;
  double Weight = 0.0;
  double VerWeight = 0.0;
  double LVerWeight, RVerWeight;
  double GWeight = 0.0;
  double TauR2L, TauL2R, SymFactor;
  int DiagIndex = 0;
  int NextLevel = Level + 1;
  int NextDiagNum = 0;

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
      VerLDiTran = InR * -1.0;

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
      int LIndex = NextDiagNum;
      Ver4Loop(VerLInL, VerLInR, VerLDiTran, loop, LTauIndex, LoopIndex + 1,
               0, // calculate u, s, t sub-ver-diagram
               Type, NextLevel, NextDiagNum);
      int LDiagNum = NextDiagNum - LIndex;

      // right vertex
      int RIndex = NextDiagNum;
      Ver4Loop(VerRInL, VerRInR, VerRDiTran, LoopNum - 1 - loop, RTauIndex,
               LoopIndex + 1 + loop,
               0, // calculate u, s, t sub-ver-diagram
               Type, NextLevel, NextDiagNum);
      int RDiagNum = NextDiagNum - RIndex;

      for (int left = LIndex; left < LIndex + LDiagNum; left++) {
        for (int right = RIndex; right < RIndex + RDiagNum; right++) {

          if ((Channel == -1 || Channel == 0) && chan == 0) {
            _ExtTau[Level][DiagNum][INL] = _ExtTau[NextLevel][left][INL];
            _ExtTau[Level][DiagNum][OUTL] = _ExtTau[NextLevel][left][OUTL];
            _ExtTau[Level][DiagNum][INR] = _ExtTau[NextLevel][right][INR];
            _ExtTau[Level][DiagNum][OUTR] = _ExtTau[NextLevel][right][OUTR];

            TauR2L =
                _ExtTau[NextLevel][left][INR] - _ExtTau[NextLevel][right][OUTL];
            TauL2R =
                _ExtTau[NextLevel][right][INL] - _ExtTau[NextLevel][left][OUTR];

          } else if ((Channel == -1 || Channel == 1) && chan == 1) {
            _ExtTau[Level][DiagNum][INL] = _ExtTau[NextLevel][left][INL];
            _ExtTau[Level][DiagNum][OUTL] = _ExtTau[NextLevel][right][OUTR];
            _ExtTau[Level][DiagNum][INR] = _ExtTau[NextLevel][right][INR];
            _ExtTau[Level][DiagNum][OUTR] = _ExtTau[NextLevel][left][OUTL];

            TauR2L =
                _ExtTau[NextLevel][left][INR] - _ExtTau[NextLevel][right][OUTL];
            TauL2R =
                _ExtTau[NextLevel][right][INL] - _ExtTau[NextLevel][left][OUTR];
          } else if ((Channel == -1 || Channel == 2) && chan == 2) {
            _ExtTau[Level][DiagNum][INL] = _ExtTau[NextLevel][left][INL];
            _ExtTau[Level][DiagNum][OUTL] = _ExtTau[NextLevel][right][OUTL];
            _ExtTau[Level][DiagNum][INR] = _ExtTau[NextLevel][left][INR];
            _ExtTau[Level][DiagNum][OUTR] = _ExtTau[NextLevel][right][OUTR];

            TauR2L =
                _ExtTau[NextLevel][left][OUTR] - _ExtTau[NextLevel][right][INR];
            TauL2R =
                _ExtTau[NextLevel][right][INL] - _ExtTau[NextLevel][left][OUTL];
          }

          VerWeight = _Weight[NextLevel][left] * _Weight[NextLevel][right];

          if (Type == -1 && _Derivative[NextLevel][left] == false &&
              _Derivative[NextLevel][right] == false) {
            // if both left and right vertex has not yet derivatived, then one
            // may appy derivative on GG
            GWeight = Fermi.Green(TauL2R, Internal2, UP, 0, Var.CurrScale) *
                      Fermi.Green(TauR2L, Internal, UP, 2, Var.CurrScale);
            GWeight += Fermi.Green(TauL2R, Internal2, UP, 2, Var.CurrScale) *
                       Fermi.Green(TauR2L, Internal, UP, 0, Var.CurrScale);
            _Derivative[Level][DiagNum] = true;
            _Weight[Level][DiagNum] = GWeight * VerWeight * (-1) * SymFactor;
            Weight += _Weight[Level][DiagNum];
            DiagNum++;
          }
          GWeight = Fermi.Green(TauL2R, Internal2, UP, 0, Var.CurrScale) *
                    Fermi.Green(TauR2L, Internal, UP, 0, Var.CurrScale);
          _Weight[Level][DiagNum] = GWeight * VerWeight * (-1) * SymFactor;
          _Derivative[Level][DiagNum] =
              _Derivative[NextLevel][left] | _Derivative[NextLevel][right];
          Weight += _Weight[Level][DiagNum];
          DiagNum++;

          // cout << "Chan:" << chan << endl;
          // cout << TauR2L << " vs " << TauL2R << endl;
          // // cout << "deltaT" << Var.Tau[0] << " vs " << Var.Tau[1] << endl;

          // cout << VerWeight << ", g: " << GWeight << ", total: " << Weight
          //      << endl;
          // cout << "Order: " << Var.CurrGroup->Order << endl;
          // cout << "LoopNum: " << LoopNum << endl;
          // cout << Var.LoopMom[1].norm() << endl;
          // cout << Var.LoopMom[3].norm() << endl;
          // cout << Internal.norm() << endl;
          // cout << Internal2.norm() << endl;
          // cout << endl;
        }
        // double VerWeight = VerLWeight[DIRECT] * VerRWeight[DIRECT];
      }
    }
  }
  return Weight;
}