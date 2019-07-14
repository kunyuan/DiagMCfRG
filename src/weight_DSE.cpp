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

cmplx weight::Evaluate(group &Group) {
  if (Group.Order == 0)
    return VerQTheta.Interaction(Var.LoopMom[1], Var.LoopMom[2], Var.LoopMom[0],
                                 -1, Var.CurrScale);
  else {
    // return fRG(group & Group);
    // RG loop
  }
}

cmplx weight::fRG(int Order) {

  const momentum &DirTran = Var.LoopMom[0];
  const momentum &InL = Var.LoopMom[1];
  const momentum &InR = Var.LoopMom[2];
  if (Order == 0) {
    // normalization
    return cmplx(VerQTheta.Interaction(InL, InR, DirTran, -2, Var.CurrScale),
                 0.0);
  } else {
    int Level = 0;
    int DiagNum = 0;
    return Ver4Loop(InL, InR, DirTran, DirTran, Order, 3,
                    1,  // calculate u diagram only
                    -1, // do RG diagram
                    Level, DiagNum);
  }
}

cmplx weight::Ver4Loop(const momentum &InL, const momentum &InR,
                       const momentum &DirTran, const momentum &ExTran,
                       int Order, int LoopIndex, int Channel, int Type,
                       int &Level, int &DiagNum) {
  cmplx Weight;
  if (Order == 0) {
    Weight = Ver4Loop0(InL, InR, DirTran, ExTran, LoopIndex, Level, DiagNum);
  } else if (Order == 1)
    Weight = Ver4Loop1(InL, InR, DirTran, ExTran, Order, LoopIndex, Channel,
                       Type, Level, DiagNum);
  return Weight;
}

cmplx weight::Ver4Loop0(const momentum &InL, const momentum &InR,
                        const momentum &DirTran, const momentum &ExTran,
                        int LoopIndex, int &Level, int &DiagNum) {
  return cmplx(VerQTheta.Interaction(InL, InR, DirTran, 0, Var.CurrScale), 0.0);
  //   return cmplx(VerQTheta.Interaction(InL, InR, DirTran, 0, Var.CurrScale),
  //                0.0) -
  //          cmplx(VerQTheta.Interaction(InL, InR, ExTran, 0, Var.CurrScale),
  //          0.0);

  //   _Weight[Level][DiagNum][DIRECT] =
  //       cmplx(VerQTheta.Interaction(InL, InR, DirTran, 0, Var.CurrScale),
  //       0.0);
  //   _Weight[Level][DiagNum][EXCHANGE] =
  //       cmplx(VerQTheta.Interaction(InL, InR, ExTran, 0, Var.CurrScale),
  //       0.0);
  //   DiagNum += 1;
  //   return _Weight[Level][DiagNum][DIRECT] -
  //   _Weight[Level][DiagNum][EXCHANGE];
}

cmplx weight::Ver4Loop1(const momentum &InL, const momentum &InR,
                        const momentum &DirTran, const momentum &ExTran,
                        int Order, int LoopIndex, int Channel, int Type,
                        int &Level, int &DiagNum) {

  momentum Internal = Var.LoopMom[LoopIndex];
  cmplx Weight = (0.0, 0.0);
  cmplx VerWeight = (0.0, 0.0);
  cmplx LVerWeight = (0.0, 0.0);
  cmplx RVerWeight = (0.0, 0.0);
  cmplx GWeight = (0.0, 0.0);
  int DiagIndex = 0;

  // u diagram
  if (Channel == 0 || Channel == 1) {
    momentum Internal2 = Internal + DirTran;
    const momentum &VerLExTran = InL - Internal - DirTran;

    const momentum &VerRInL = Internal2;
    const momentum &VerRExTran = Internal - InR;

    for (int subOrder = 0; subOrder < Order; subOrder++) {
      //====================  DIRECT  Diagram =============================
      // left vertex
      int LIndex = DiagNum;
      LVerWeight +=
          Ver4Loop(InL, Internal, DirTran, VerLExTran, subOrder, LoopIndex + 1,
                   0, // calculate u, s, t sub-ver-diagram
                   0, // type 0
                   Level, DiagNum);
      int LDiagNum = DiagNum - LIndex;

      // right vertex
      int RIndex = DiagNum;
      RVerWeight += Ver4Loop(VerRInL, InR, DirTran, VerRExTran,
                             Order - 1 - subOrder, LoopIndex + 1 + subOrder,
                             0, // calculate u, s, t sub-ver-diagram
                             0, // type 0
                             Level, DiagNum);
      int RDiagNum = DiagNum - RIndex;
    }
    if (Type == -1) {
      GWeight = Fermi.GreenFreq(Internal2, 0, Var.CurrScale) *
                Fermi.GreenFreq(Internal, 2, Var.CurrScale);
      GWeight += Fermi.GreenFreq(Internal2, 2, Var.CurrScale) *
                 Fermi.GreenFreq(Internal, 0, Var.CurrScale);
    }

    Weight += GWeight * LVerWeight * RVerWeight * pow(-1, Order);

    // cout << "Order" << Order << endl;
    // cout << LVerWeight << " ver " << RVerWeight << endl;
    // cout << VerWeight << ", g: " << GWeight << ", total: " << Weight << endl;
    // cout << "inloop2: " << Internal2.norm() << " freq: " << Internal2[D]
    //      << endl;
    // cout << "inloop1: " << Internal.norm() << " freq: " << Internal[D] <<
    // endl;
  }
  //   Weight *= pow(-1, LoopNum);
  return Weight;
}