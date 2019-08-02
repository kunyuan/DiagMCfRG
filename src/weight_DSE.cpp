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

#define SETTAU(Order, Level, Index, tInL, tOUTL, tINR, tOUTR)                  \
  if (_ExtTau[Order][Level][Index][INL] == -1) {                               \
    _ExtTau[Order][Level][Index][INL] = tInL;                                  \
    _ExtTau[Order][Level][Index][OUTL] = tOUTL;                                \
    _ExtTau[Order][Level][Index][INR] = tINR;                                  \
    _ExtTau[Order][Level][Index][OUTR] = tOUTR;                                \
  }

double weight::Evaluate(group &Group) {
  if (Group.Order == 0)
    VerQTheta.Interaction(Var.LoopMom[1], Var.LoopMom[2], Var.LoopMom[0], -1,
                          Var.CurrScale);
  else {
    // return fRG(group & Group);
    // RG loop
  }
}

double weight::fRG(int Order, int ID) {

  const momentum &DirTran = Var.LoopMom[0];
  const momentum &InL = Var.LoopMom[1];
  const momentum &InR = Var.LoopMom[2];
  if (Order == 0) {
    // normalization
    _DiagIndex[0] = 1;
    _Weight[0][0][0] = VerQTheta.Interaction(InL, InR, DirTran, 0.0, -2);
    // _ExtTau[0][0][INL]=0.0;
    return _Weight[0][0][0];
  } else {
    _GlobalOrder = Order;

    for (int i = 0; i < Order + 1; i++)
      _DiagIndex[i] = 0;

    int Level = 0;
    if (Order == 1) {
      Vertex4(InL, InR, DirTran, Order, 0, 3, Level,
              T,  // t diagram only
              -1, // normal diagram
              -1  // left vertex order
      );
      // if (ID == 0)
      //   return _Weight[0][0][0];
      // Weight = _Weight[Level][0];
    } else if (Order == 2) {
      Vertex4(InL, InR, DirTran, Order, 0, 3, Level,
              T,  // t diagram only
              -1, // normal diagram
                  // RIGHT, // normal diagram
              -1  // left vertex order
                  // 1 // left vertex order
      );
    } else if (Order == 3) {
      Vertex4(InL, InR, DirTran, Order, 0, 3, Level,
              T,  // t diagram only
              -1, // normal diagram
              // 1, // diff diagram
              -1 // left vertex order
      );
    }
    double Weight = 0.0;
    int count = 0;
    // double inL = _ExtTau[Level][0][INL];
    for (int diag = 0; diag < _DiagIndex[Level]; diag++) {
      Weight += _Weight[Level][diag][0];
      count++;
      // if (abs(_ExtTau[Level][diag][INL] - inL) > 1.0e-10) {
      //   ABORT("left time different! DiagNum: " << diag << " with " << inL
      //                                          << " vs "
      //                                          << _ExtTau[Level][diag][INL]);
      // }
    }
    _DiagNum = count;
    // if (LoopNum == 3 || LoopNum == 2)
    //   cout << "order: " << LoopNum << ", RG number: " << count
    //        << ", diag num: " << DiagIndex << endl;

    // cout << _Weight[Level][0] << ": " << _Weight[Level][1] << ": "
    //      << _Weight[Level][2] << endl;

    // if (Order == 2)
    //   // cout << Weight << endl;
    // cout << count << endl;
    return Weight / pow(2.0 * PI, 2 * Order);
    // return Weight;
  }
}

int weight::Vertex4(
    const momentum &InL, const momentum &InR, const momentum &DirTran,
    int Order, int TauIndex, int LoopIndex, int Level,
    bool *Channel, // three flags, calculate t, u, s or not
    int VerType,   // -1: normal, 0: left(to project), 1: right(to diff)
    int LVerOrder  // order of left vertex
) {
  if (Order == 0) {
    return Ver4Loop0(InL, InR, DirTran, TauIndex, LoopIndex, Level);
  } else if (Order >= 1) {

    // if (LoopNum == 2 && Level == 1)
    //   cout << "1 order: " << LoopNum << "diag:" << DiagIndex << endl;

    Bubble(InL, InR, DirTran, Order, TauIndex, LoopIndex, Level, Channel,
           VerType,   // VerType
           LVerOrder, // no projection
                      //  0,       // no projection
           false      // not penguin diagram
    );

    // if (LoopNum == 2 && Level == 0)
    //   cout << "2 order: " << LoopNum << "diag:" << DiagIndex << endl;

    if (Order >= 2) {
      Bubble(InL, InR, DirTran, Order, TauIndex, LoopIndex, Level, Channel,
             VerType,   // VerType
             LVerOrder, // no projection
             true       // penguin diagram
      );
    }
    // if (LoopNum == 2 && Level == 0)
    //   cout << "2 order: " << LoopNum << "diag:" << DiagIndex << endl;
    // for normal vertex or projected vertex, just return
    // penguin diagram
    return _DiagIndex[Level];
  }
}

int weight::Bubble(
    const momentum &InL, const momentum &InR, const momentum &DirTran,
    int Order, int TauIndex, int LoopIndex, int Level,
    bool *Channel, // three flags, calculate t, u, s or not
    int VerType,   // -1: normal, 0: left(to project), 1: right(to diff)
    int LVerOrder, // order of left vertex
    bool IsPenguin) {

  // calculate renormalized diagrams
  for (int OL = 0; OL < Order; OL++) {
    if (LVerOrder >= 0 && OL != LVerOrder)
      continue;
    if (IsPenguin && OL < 1)
      continue;

    if (VerType == -1 || VerType == 1) {
      // for normal vertex or projected vertex, just return
      OneLoop(InL, InR, DirTran, Order, OL, TauIndex, LoopIndex, Level, Channel,
              false, // do not project
              IsPenguin);
      // if (Order == 2 && Level == 0) {
      //   cout << VerType << ", index=" << _DiagIndex[Level]
      //        << ", level=" << Level << " OL:" << OL << endl;
      // }
    }

    /////// comment this block to calculate original diagrams  //////
    if (VerType == 0 || VerType == 1) {
      // do projection
      OneLoop(InL, InR, DirTran, Order, OL, TauIndex, LoopIndex, Level, Channel,
              true, // do projection
              IsPenguin);
      // if (Order == 2 && Level == 0) {
      //   cout << VerType << ", index=" << _DiagIndex[Level]
      //        << ", level=" << Level << " OL:" << OL << endl;
      // }
    }
  }
  return _DiagIndex[Level];
}

int weight::Ver4Loop0(const momentum &InL, const momentum &InR,
                      const momentum &DirTran, int TauIndex, int LoopIndex,
                      int Level, int Type) {
  //   return VerQTheta.Interaction(InL, InR, DirTran, 0, Var.CurrScale) -
  //   VerQTheta.Interaction(InL, InR, ExTran, 0, Var.CurrScale);
  double Tau = Var.Tau[TauIndex] - Var.Tau[TauIndex + 1];
  int &Index = _DiagIndex[Level];
  ////////////// bare interaction ///////////
  double DiWeight = VerQTheta.Interaction(InL, InR, DirTran, Tau, 0);
  double ExWeight =
      VerQTheta.Interaction(InL, InR, InR + DirTran - InL, Tau, 0);
  SETTAU(_GlobalOrder, Level, Index, TauIndex, TauIndex, TauIndex, TauIndex);
  _Weight[Level][Index][0] = DiWeight - ExWeight;
  // _Weight[Level][DiagIndex][0] = DiWeight;
  Index += 1;

  if (Type != -2) {
    //////////// dressed interaction ///////////
    DiWeight = VerQTheta.Interaction(InL, InR, DirTran, Tau, 1);
    SETTAU(_GlobalOrder, Level, Index, TauIndex, TauIndex, TauIndex + 1,
           TauIndex + 1);
    _Weight[Level][Index][0] = DiWeight;
    Index += 1;

    ExWeight = VerQTheta.Interaction(InL, InR, InR + DirTran - InL, Tau, 1);
    SETTAU(_GlobalOrder, Level, Index, TauIndex, TauIndex + 1, TauIndex + 1,
           TauIndex);
    _Weight[Level][Index][0] = -ExWeight;
    Index += 1;
  }

  return Index;
}

int weight::OneLoop(const momentum &InL, const momentum &InR,
                    const momentum &DirTran, int Order, int LVerOrder,
                    int TauIndex, int LoopIndex, int Level,
                    bool *Channel, // three flags, calculate t, u, s or not
                    bool IsProjected, bool IsPenguin) {

  double VerWeight;
  double GR2L, GL2R;
  double SymFactor = 1.0;
  int nLevel = Level + 1;
  int &Index = _DiagIndex[Level];

  if (Order < 1)
    return Index;

  if (IsPenguin && (LVerOrder < 1 || Order < 2))
    return Index;

  // do projection
  double ProjSign = 1.0;
  if (IsProjected)
    ProjSign = -1.0; // projection always comes with a minus sign

  momentum Internal = Var.LoopMom[LoopIndex + LVerOrder];
  momentum Internal2, VerLInL, VerLInR, VerLDiTran, VerRInL, VerRInR,
      VerRDiTran;

  int LTauIndex = TauIndex;
  int RTauIndex = TauIndex + (LVerOrder + 1) * 2;

  for (int chan = 0; chan < 3; chan++) {
    if (!Channel[chan])
      continue;
    if (chan == 0) {
      // t diagram
      if (IsProjected) {
        VerLInL = InL * (Para.Kf / InL.norm());
        VerRInR = InR * (Para.Kf / InR.norm());
      } else {
        VerLInL = InL;
        VerRInR = InR;
      }
      Internal2 = Internal - DirTran;
      VerLInR = Internal2;
      VerLDiTran = DirTran;

      VerRInL = Internal;
      VerRDiTran = DirTran;

      SymFactor = -1.0;
    } else if (chan == 1) {
      // u diagram
      if (IsProjected) {
        VerLInL = InL * (Para.Kf / InL.norm());
        VerRInR = InR * (Para.Kf / InR.norm());
      } else {
        VerLInL = InL;
        VerRInR = InR;
      }
      Internal2 = Internal + DirTran - InL + InR;
      // after the projection, the exchange transfer momentum
      // should remain the same
      VerLInR = Internal2;
      VerLDiTran = Internal - Internal2;
      // after the projection, the direct transfer momentum on left and right
      // vertex also remain the same !!!

      VerRInL = Internal;
      VerRDiTran = VerLDiTran;
      SymFactor = 1.0;
    } else if (chan == 2) {
      // projection is non-zero only for t and u channel
      if (IsProjected)
        continue;
      // if (IsProjected) {
      //   VerLInL = InL * (Para.Kf / InL.norm());
      //   VerRInR = InR * (Para.Kf / InR.norm());
      // } else {
      //   VerLInL = InL;
      //   VerRInR = InR;
      // }

      // s diagram
      Internal2 = InL + InR - Internal;
      VerLInL = InL;
      VerLInR = InR;
      VerLDiTran = VerLInL - Internal2;

      VerRInL = Internal2;
      VerRInR = Internal;
      VerRDiTran = DirTran + InR - Internal;
      SymFactor = 0.5;
    }

    //====================  DIRECT  Diagram =============================
    // left vertex
    bool *nChannel = ALL;
    int LVerType = LEFT;
    int LIndex = _DiagIndex[nLevel];

    if (IsPenguin) {
      LVerType = -1; // normal left vertex;
      if (chan == 0 || chan == 1)
        nChannel = US;
      else
        nChannel = UT;
      Bubble(VerLInL, VerLInR, VerLDiTran, LVerOrder, LTauIndex, LoopIndex,
             nLevel, nChannel,
             LVerType, // VerType
             -1,       // LVerOrder
             false     // not penguin
      );
    } else {
      // nChannel = T;
      // if (Level == 0)
      //   LVerType = -2;
      Vertex4(VerLInL, VerLInR, VerLDiTran, LVerOrder, LTauIndex, LoopIndex,
              nLevel, nChannel,
              LVerType, // VerType
              -1        // LVerOrder
      );
    }
    int LDiagIndex = _DiagIndex[nLevel];

    // nChannel = T;
    // right vertex
    int RIndex = LDiagIndex;
    Vertex4(VerRInL, VerRInR, VerRDiTran, Order - 1 - LVerOrder, RTauIndex,
            LoopIndex + 1 + LVerOrder, nLevel, ALL,
            // nChannel,
            RIGHT, // VerType
            // -1, // VerType
            -1 // LVerOrder
    );

    // if (LVerLoopNum == 0 && LoopNum == 2) {
    //   nChannel = U;
    //   // right vertex
    //   nDiagIndex =
    //       Vertex4(VerRInL, VerRInR, VerRDiTran, LoopNum - 1 - LVerLoopNum,
    //               RTauIndex, LoopIndex + 1 + LVerLoopNum, nDiagIndex, nLevel,
    //               // ALL,
    //               nChannel,
    //               RIGHT, // VerType
    //               -1     // LVerOrder
    //       );
    // }
    int RDiagIndex = _DiagIndex[nLevel];

    ////////////// construct  G table  /////////////////////////////////////
    // for (int tL = LTauIndex; tL < RTauIndex; tL++)
    //   for (int tR = RTauIndex; tR < TauIndex + (Order + 1) * 2; tR++) {
    //     _GL2R[tL][tR] = Fermi.Green(Var.Tau[tR] - Var.Tau[tL], Internal, UP,
    //     0,
    //                                 Var.CurrScale);
    //     if (chan == 2)
    //       _GR2L[tL][tR] = Fermi.Green(Var.Tau[tR] - Var.Tau[tL], Internal2,
    //       UP,
    //                                   0, Var.CurrScale);
    //     else
    //       _GR2L[tR][tL] = Fermi.Green(Var.Tau[tL] - Var.Tau[tR], Internal2,
    //       UP,
    //                                   0, Var.CurrScale);
    //     // cout << tR << "-" << tL << ", " << _GR2L[tR][tL] << endl;
    //   }

    if (chan == 2) {
      for (int tL = LTauIndex; tL < RTauIndex; tL++) {
        _GR2L[tL][RTauIndex] = Fermi.Green(Var.Tau[RTauIndex] - Var.Tau[tL],
                                           Internal2, UP, 0, Var.CurrScale);
        for (int tR = RTauIndex; tR < TauIndex + (Order + 1) * 2; tR++)
          _GL2R[tL][tR] = Fermi.Green(Var.Tau[tR] - Var.Tau[tL], Internal, UP,
                                      0, Var.CurrScale);
      }
    } else {
      for (int tL = LTauIndex; tL < RTauIndex; tL++) {
        _GL2R[tL][RTauIndex] = Fermi.Green(Var.Tau[RTauIndex] - Var.Tau[tL],
                                           Internal, UP, 0, Var.CurrScale);
        for (int tR = RTauIndex; tR < TauIndex + (Order + 1) * 2; tR++)
          _GR2L[tR][tL] = Fermi.Green(Var.Tau[tL] - Var.Tau[tR], Internal2, UP,
                                      0, Var.CurrScale);
      }
    }

    int *_ExtTauC = _ExtTau[_GlobalOrder][Level][Index];
    for (int l = LIndex; l < LDiagIndex; l++) {
      int *_ExtTauL = _ExtTau[_GlobalOrder][nLevel][l];
      for (int r = RIndex; r < RDiagIndex; r++) {
        int *_ExtTauR = _ExtTau[_GlobalOrder][nLevel][r];

        if (chan == 0) {
          SETTAU(_GlobalOrder, Level, Index, _ExtTauL[INL], _ExtTauL[OUTL],
                 _ExtTauR[INR], _ExtTauR[OUTR]);

          GR2L = _GR2L[_ExtTauR[OUTL]][_ExtTauL[INR]];
          GL2R = _GL2R[_ExtTauL[OUTR]][_ExtTauR[INL]];
          // cout << "0me: " << GR2L << ", " << GL2R << endl;
        } else if (chan == 1) {
          SETTAU(_GlobalOrder, Level, Index, _ExtTauL[INL], _ExtTauR[OUTR],
                 _ExtTauR[INR], _ExtTauL[OUTL]);

          GR2L = _GR2L[_ExtTauR[OUTL]][_ExtTauL[INR]];
          GL2R = _GL2R[_ExtTauL[OUTR]][_ExtTauR[INL]];
          // cout << "1me: " << GR2L << ", " << GL2R << endl;
        } else if (chan == 2) {
          SETTAU(_GlobalOrder, Level, Index, _ExtTauL[INL], _ExtTauR[OUTL],
                 _ExtTauL[INR], _ExtTauR[OUTR]);

          GR2L = _GR2L[_ExtTauL[OUTL]][_ExtTauR[INL]];
          GL2R = _GL2R[_ExtTauL[OUTR]][_ExtTauR[INR]];
          // cout << "2me: " << GR2L << ", " << GL2R << endl;
        }

        if (IsProjected)
          if (chan == 0) {
            SETTAU(_GlobalOrder, Level, Index, _ExtTauC[INL], _ExtTauC[INL],
                   _ExtTauC[INR], _ExtTauC[INR]);
          } else if (chan == 1) {
            SETTAU(_GlobalOrder, Level, Index, _ExtTauC[INL], _ExtTauC[INR],
                   _ExtTauC[INR], _ExtTauC[INL]);
          }

        _Weight[Level][Index][0] = _Weight[nLevel][l][0] *
                                   _Weight[nLevel][r][0] * GL2R * GR2L *
                                   SymFactor * ProjSign;

        // if (_GlobalOrder == 2 && l == LIndex && r == RIndex) {
        //   cout << "Lindex: " << LIndex << " Weight:" <<
        //   _Weight[Level][Index][0]
        //        << endl;
        //   cout << Var.Tau[_ExtTauL[OUTR]] << ", " << Var.Tau[_ExtTauR[INL]]
        //        << ", Mom=" << Internal.norm() << ": GL2R=" << GL2R << endl;
        //   // cout << Fermi.Green(Var.Tau[_ExtTauR[INL]] -
        //   // Var.Tau[_ExtTauL[OUTR]],
        //   //                     Internal, UP, 0, Var.CurrScale)
        //   //      << endl;
        //   cout << _ExtTauR[OUTL] << ", " << _ExtTauL[INR]
        //        << ", Mom=" << Internal2.norm() << ": GR2L=" << GR2L
        //        << ", chan=" << chan << endl;
        //   cout << Fermi.Green(Var.Tau[_ExtTauL[INR]] -
        //   Var.Tau[_ExtTauR[OUTL]],
        //                       Internal2, UP, 0, Var.CurrScale)
        //        << endl;
        //   cout << _Weight[nLevel][l][0] << endl;
        //   cout << _Weight[nLevel][r][0] << endl;
        //   // cout << GL2R << " vs " << GR2L << endl;
        //   cout << "end" << endl;
        //   cout << endl;
        // }

        Index++;
      }
    }
  }
  return Index;
}