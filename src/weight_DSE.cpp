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
    return VerQTheta.Interaction(InL, InR, DirTran, 0.0, -2);
  } else {
    int Level = 0;
    int DiagIndex = 0;
    double Weight = 0.0;
    if (LoopNum == 1) {
      DiagIndex = Vertex4(InL, InR, DirTran, LoopNum, 0, 3, DiagIndex, Level,
                          T,  // t diagram only
                          -1, // normal diagram
                          -1  // left vertex order
      );
      // Weight = _Weight[Level][0];
    } else if (LoopNum == 2) {
      DiagIndex = Vertex4(InL, InR, DirTran, LoopNum, 0, 3, DiagIndex, Level,
                          T,  // t diagram only
                          -1, // normal diagram
                          // 1, // diff diagram
                          -1 // left vertex order
      );
    } else if (LoopNum == 3) {
      DiagIndex = Vertex4(InL, InR, DirTran, LoopNum, 0, 3, DiagIndex, Level,
                          T,  // t diagram only
                          -1, // normal diagram
                          // 1, // diff diagram
                          -1 // left vertex order
      );
    }
    int count = 0;
    for (int diag = 0; diag < DiagIndex; diag++) {
      Weight += _Weight[Level][diag][1];
      count++;
    }
    // if (LoopNum == 3 || LoopNum == 2)
    //   cout << "order: " << LoopNum << ", RG number: " << count
    //        << ", diag num: " << DiagIndex << endl;

    // cout << _Weight[Level][0] << ": " << _Weight[Level][1] << ": "
    //      << _Weight[Level][2] << endl;

    return Weight / pow(40.0, LoopNum);
  }
}

int weight::Vertex4(
    const momentum &InL, const momentum &InR, const momentum &DirTran,
    int LoopNum, int TauIndex, int LoopIndex, int DiagIndex, int Level,
    bool *Channel, // three flags, calculate t, u, s or not
    int VerType,   // -1: normal, 0: left(to project), 1: right(to diff)
    int LVerOrder  // order of left vertex
) {
  if (LoopNum == 0) {
    return Ver4Loop0(InL, InR, DirTran, TauIndex, LoopIndex, Level, DiagIndex);
  } else if (LoopNum >= 1) {

    // if (LoopNum == 2 && Level == 1)
    //   cout << "1 order: " << LoopNum << "diag:" << DiagIndex << endl;

    DiagIndex = Bubble(InL, InR, DirTran, LoopNum, TauIndex, LoopIndex,
                       DiagIndex, Level, Channel,
                       VerType,   // VerType
                       LVerOrder, // no projection
                       false      // not penguin diagram
    );

    // if (LoopNum == 2 && Level == 1)
    //   cout << "2 order: " << LoopNum << "diag:" << DiagIndex << endl;

    if (LoopNum >= 2) {
      DiagIndex = Bubble(InL, InR, DirTran, LoopNum, TauIndex, LoopIndex,
                         DiagIndex, Level, Channel,
                         VerType,   // VerType
                         LVerOrder, // no projection
                         true       // penguin diagram
      );
    }
    // if (LoopNum == 2 && Level == 1)
    //   cout << "3 order: " << LoopNum << "diag:" << DiagIndex << endl;
    // for normal vertex or projected vertex, just return
    // penguin diagram
    return DiagIndex;
  }
}

int weight::Bubble(
    const momentum &InL, const momentum &InR, const momentum &DirTran,
    int LoopNum, int TauIndex, int LoopIndex, int DiagIndex, int Level,
    bool *Channel, // three flags, calculate t, u, s or not
    int VerType,   // -1: normal, 0: left(to project), 1: right(to diff)
    int LVerOrder, // order of left vertex
    bool IsPenguin) {

  // calculate normal diagrams
  // for (int OL = 0; OL < LoopNum; OL++) {
  //   if (IsPenguin && OL < 1)
  //     continue;
  //   if (!IsPenguin && OL > 0)
  //     continue;
  //   DiagIndex = OneLoop(InL, InR, DirTran, LoopNum, OL, TauIndex, LoopIndex,
  //                       DiagIndex, Level, Channel,
  //                       false, // do not project
  //                       IsPenguin);
  // }

  // calculate renormalized diagrams
  for (int OL = 0; OL < LoopNum; OL++) {
    if (LVerOrder >= 0 && OL != LVerOrder)
      continue;
    if (IsPenguin && OL < 1)
      continue;

    if (VerType == -1 || VerType == 1) {
      // for normal vertex or projected vertex, just return
      DiagIndex = OneLoop(InL, InR, DirTran, LoopNum, OL, TauIndex, LoopIndex,
                          DiagIndex, Level, Channel,
                          false, // do not project
                          IsPenguin);
      // if (LoopNum == 2) {
      //   cout << VerType << ", index=" << DiagIndex << ", level=" << Level
      //        << " OL:" << OL << endl;
      // }
    }
    if (VerType == 0 || VerType == 1) {
      // do projection
      DiagIndex = OneLoop(InL, InR, DirTran, LoopNum, OL, TauIndex, LoopIndex,
                          DiagIndex, Level, Channel,
                          true, // do projection
                          IsPenguin);
    }
  }
  return DiagIndex;
}

int weight::Ver4Loop0(const momentum &InL, const momentum &InR,
                      const momentum &DirTran, int TauIndex, int LoopIndex,
                      int Level, int DiagIndex) {
  //   return VerQTheta.Interaction(InL, InR, DirTran, 0, Var.CurrScale) -
  //   VerQTheta.Interaction(InL, InR, ExTran, 0, Var.CurrScale);
  double Tau = Var.Tau[TauIndex] - Var.Tau[TauIndex + 1];
  double DiWeight = -VerQTheta.Interaction(InL, InR, DirTran, abs(Tau), 0);
  double ExWeight =
      -VerQTheta.Interaction(InL, InR, InR + DirTran - InL, abs(Tau), 0);
  _ExtTau[Level][DiagIndex][INL] = Var.Tau[TauIndex];
  _ExtTau[Level][DiagIndex][OUTL] = Var.Tau[TauIndex];
  _ExtTau[Level][DiagIndex][INR] = Var.Tau[TauIndex + 1];
  _ExtTau[Level][DiagIndex][OUTR] = Var.Tau[TauIndex + 1];
  _Weight[Level][DiagIndex][0] = DiWeight - ExWeight;
  // _Weight[Level][DiagIndex][0] = DiWeight;
  _Weight[Level][DiagIndex][1] = 0.0;

  // ASSERT_ALLWAYS(abs(DirTran.norm() - Var.LoopMom[0].norm()) < 1.0e-5,
  //                "Ext Mom wrong!");

  DiagIndex += 1;
  return DiagIndex;
}

int weight::OneLoop(const momentum &InL, const momentum &InR,
                    const momentum &DirTran, int LoopNum, int LVerLoopNum,
                    int TauIndex, int LoopIndex, int DiagIndex, int Level,
                    bool *Channel, // three flags, calculate t, u, s or not
                    bool IsProjected, bool IsPenguin) {

  double VerWeight;
  double TauR2L, TauL2R;
  double SymFactor = 1.0;
  int nLevel = Level + 1;
  int nDiagIndex = 0;

  if (LoopNum < 1)
    return DiagIndex;

  if (IsPenguin && (LVerLoopNum < 1 || LoopNum < 2))
    return DiagIndex;

  // do projection
  double ProjSign = 1.0;
  if (IsProjected)
    ProjSign = -1.0; // projection always comes with a minus sign

  momentum Internal = Var.LoopMom[LoopIndex + LVerLoopNum];
  momentum Internal2, VerLInL, VerLInR, VerLDiTran, VerRInL, VerRInR,
      VerRDiTran;

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

    int LTauIndex = TauIndex;
    int RTauIndex = TauIndex + (LVerLoopNum + 1);

    //====================  DIRECT  Diagram =============================
    // left vertex
    bool *nChannel = ALL;
    int LVerType = LEFT;
    int LIndex = nDiagIndex;
    if (IsPenguin) {
      LVerType = -1; // normal left vertex;
      if (chan == 0 || chan == 1)
        nChannel = US;
      else
        nChannel = UT;
      nDiagIndex = Bubble(VerLInL, VerLInR, VerLDiTran, LVerLoopNum, LTauIndex,
                          LoopIndex, nDiagIndex, nLevel, nChannel,
                          LVerType, // VerType
                          -1,       // LVerOrder
                          false     // not penguin
      );
    } else {
      nChannel = T;
      nDiagIndex = Vertex4(VerLInL, VerLInR, VerLDiTran, LVerLoopNum, LTauIndex,
                           LoopIndex, nDiagIndex, nLevel, nChannel,
                           LVerType, // VerType
                           -1        // LVerOrder
      );
    }
    int LDiagIndex = nDiagIndex;

    // nChannel = T;
    // right vertex
    int RIndex = LDiagIndex;
    // if (LoopNum == 2 && LVerLoopNum == 0) {
    //   int a = 2;
    // }
    nDiagIndex =
        Vertex4(VerRInL, VerRInR, VerRDiTran, LoopNum - 1 - LVerLoopNum,
                RTauIndex, LoopIndex + 1 + LVerLoopNum, nDiagIndex, nLevel, ALL,
                // nChannel,
                RIGHT, // VerType
                -1     // LVerOrder
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
    int RDiagIndex = nDiagIndex;

    for (int left = LIndex; left < LDiagIndex; left++) {
      for (int right = RIndex; right < RDiagIndex; right++) {
        if (chan == 0) {
          _ExtTau[Level][DiagIndex][INL] = _ExtTau[nLevel][left][INL];
          _ExtTau[Level][DiagIndex][OUTL] = _ExtTau[nLevel][left][OUTL];
          _ExtTau[Level][DiagIndex][INR] = _ExtTau[nLevel][right][INR];
          _ExtTau[Level][DiagIndex][OUTR] = _ExtTau[nLevel][right][OUTR];
        } else if (chan == 1) {
          _ExtTau[Level][DiagIndex][INL] = _ExtTau[nLevel][left][INL];
          _ExtTau[Level][DiagIndex][OUTL] = _ExtTau[nLevel][right][OUTR];
          _ExtTau[Level][DiagIndex][INR] = _ExtTau[nLevel][right][INR];
          _ExtTau[Level][DiagIndex][OUTR] = _ExtTau[nLevel][left][OUTL];
        } else if (chan == 2) {
          _ExtTau[Level][DiagIndex][INL] = _ExtTau[nLevel][left][INL];
          _ExtTau[Level][DiagIndex][OUTL] = _ExtTau[nLevel][right][OUTL];
          _ExtTau[Level][DiagIndex][INR] = _ExtTau[nLevel][left][INR];
          _ExtTau[Level][DiagIndex][OUTR] = _ExtTau[nLevel][right][OUTR];
        }
        if (IsProjected == false) {
        } else {
          // double avg = _ExtTau[nLevel][left][INL];
          // _ExtTau[Level][DiagIndex][INL] = avg;
          // _ExtTau[Level][DiagIndex][OUTL] = avg;
          // _ExtTau[Level][DiagIndex][INR] = avg;
          // _ExtTau[Level][DiagIndex][OUTR] = avg;

          // double avg = (_ExtTau[Level][DiagIndex][INL] +
          //               _ExtTau[Level][DiagIndex][OUTL] +
          //               _ExtTau[Level][DiagIndex][INR] +
          //               _ExtTau[Level][DiagIndex][OUTR]) /
          //              4.0;
          // _ExtTau[Level][DiagIndex][INL] = avg;
          // _ExtTau[Level][DiagIndex][OUTL] = avg;
          // _ExtTau[Level][DiagIndex][INR] = avg;
          // _ExtTau[Level][DiagIndex][OUTR] = avg;

          if (chan == 0) {
            double avg = _ExtTau[Level][DiagIndex][INL];
            _ExtTau[Level][DiagIndex][INL] = avg;
            _ExtTau[Level][DiagIndex][OUTL] = avg;
            avg = _ExtTau[Level][DiagIndex][INR];
            _ExtTau[Level][DiagIndex][INR] = avg;
            _ExtTau[Level][DiagIndex][OUTR] = avg;
          } else if (chan == 1) {
            double avg = _ExtTau[Level][DiagIndex][INL];
            _ExtTau[Level][DiagIndex][INL] = avg;
            _ExtTau[Level][DiagIndex][OUTR] = avg;
            avg = _ExtTau[Level][DiagIndex][INR];
            _ExtTau[Level][DiagIndex][INR] = avg;
            _ExtTau[Level][DiagIndex][OUTL] = avg;
          }
        }

        if (chan == 2) {
          TauR2L = _ExtTau[nLevel][right][INL] - _ExtTau[nLevel][left][OUTL];
          TauL2R = _ExtTau[nLevel][right][INR] - _ExtTau[nLevel][left][OUTR];
        } else {
          TauR2L = _ExtTau[nLevel][left][INR] - _ExtTau[nLevel][right][OUTL];
          TauL2R = _ExtTau[nLevel][right][INL] - _ExtTau[nLevel][left][OUTR];
        }

        _Weight[Level][DiagIndex][0] = 0.0;
        _Weight[Level][DiagIndex][1] = 0.0;

        double GL2R = Fermi.Green(TauL2R, Internal, UP, 0, Var.CurrScale);
        double GR2L = Fermi.Green(TauR2L, Internal2, UP, 0, Var.CurrScale);

        VerWeight = _Weight[nLevel][left][0] * _Weight[nLevel][right][0];
        _Weight[Level][DiagIndex][0] +=
            GL2R * GR2L * VerWeight * SymFactor * ProjSign;

        if (IsProjected == false) {
          // if projected, then derivative must be zero!!!
          // take derivative
          double GL2RDeri = Fermi.Green(TauL2R, Internal, UP, 2, Var.CurrScale);
          double GR2LDeri =
              Fermi.Green(TauR2L, Internal2, UP, 2, Var.CurrScale);
          // if both left and right vertex do not contain
          // derivative, then GG can contain derivative
          VerWeight = _Weight[nLevel][left][0] * _Weight[nLevel][right][0];
          _Weight[Level][DiagIndex][1] += (GL2R * GR2LDeri + GL2RDeri * GR2L) *
                                          VerWeight * SymFactor * ProjSign;

          // one of the vertex function contain derivative
          VerWeight = _Weight[nLevel][left][0] * _Weight[nLevel][right][1];
          VerWeight += _Weight[nLevel][left][1] * _Weight[nLevel][right][0];
          _Weight[Level][DiagIndex][1] +=
              GL2R * GR2L * VerWeight * SymFactor * ProjSign;
        }
        DiagIndex++;
      }
    }
  }
  return DiagIndex;
}