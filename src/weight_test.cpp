#define FMT_HEADER_ONLY
#include "utility/fmt/format.h"
#include "utility/utility.h"
#include "weight.h"
#include <cmath>
#include <iostream>

using namespace diag;
using namespace std;

string weight::DebugInfo(group &Group) {
  string msg;
  msg = string(80, '=') + "\n";
  msg += fmt::format("\nMC Counter: {} \n", Para.Counter);
  msg += "Current Group Info:\n ";
  msg += ToString(Group);

  // msg += string(80, '=') + "\n";
  // msg += "GWeight: \n";
  // for (int i = 0; i < Group.GNum; i++)
  //   msg += ToString(Group.Diag[0].G[i]->Weight) + "; ";
  // msg += "\n";

  // msg += "NewGWeight: \n";
  // for (int i = 0; i < Group.GNum; i++)
  //   msg += ToString(Group.Diag[0].G[i]->NewWeight) + "; ";
  // msg += "\n";

  // msg += "TauBasis: \n";
  // for (int i = 0; i < Var.CurrGroup->GNum; i++) {
  //   msg += ToString(Var.CurrGroup->Diag[0].G[i]->TauBasis[0]) + "; ";
  // }
  // msg += "\n";
  // for (int i = 0; i < Var.CurrGroup->GNum; i++) {
  //   msg += ToString(Var.CurrGroup->Diag[0].G[i]->TauBasis[1]) + "; ";
  // }
  // msg += "\n";

  // msg += "LoopBasis: \n";
  // for (int i = 0; i < Var.CurrGroup->GNum; i++) {
  //   msg += ToString(Var.CurrGroup->Diag[0].G[i]->LoopBasis[0]) + "; ";
  // }
  // msg += "\n";
  // for (int i = 0; i < Var.CurrGroup->GNum; i++) {
  //   msg += ToString(Var.CurrGroup->Diag[0].G[i]->LoopBasis[1]) + "; ";
  // }
  // msg += "\n";

  // double Tau = Var.Tau[1] - Var.Tau[0];
  // momentum G1, G2;
  // for (int i = 0; i < D; i++) {
  //   G1[i] = Var.LoopMom[0][i] + Var.LoopMom[1][i];
  //   G2[i] = Var.LoopMom[1][i];
  // }
  // msg += "Expected G: \n" + ToString(Green(Tau, G1, UP, 0, true)) + "; " +
  //        ToString(Green(-Tau, G2, UP, 0, true)) + "; ";
  // msg += "\n";

  // msg += "Ver4Weight: \n";
  // for (int i = 0; i < Group.Ver4Num; i++)
  //   msg += fmt::format("({:.3f}, {:.3f}); ",
  //   Group.Diag[0].Ver4[i]->Weight[0],
  //                      Group.Diag[0].Ver4[i]->Weight[1]);
  // msg += "\n";

  // msg += "NewVer4Weight: \n";
  // for (int i = 0; i < Group.Ver4Num; i++)
  //   msg +=
  //       fmt::format("({:.3f}, {:.3f}); ",
  //       Group.Diag[0].Ver4[i]->NewWeight[0],
  //                   Group.Diag[0].Ver4[i]->NewWeight[1]);
  // msg += "\n";

  msg += string(80, '=') + "\n";
  msg += "LoopMom: \n";
  for (int i = 0; i < Group.LoopNum; i++)
    msg +=
        ToString(Var.LoopMom[i]) + "=" + ToString(Var.LoopMom[i].norm()) + "\n";

  msg += string(80, '=') + "\n";
  msg += "Tau: \n";
  for (int i = 0; i < Group.TauNum; i++)
    msg += ToString(Var.Tau[i]) + ", ";

  msg += "\n";
  return msg;
}

string weight::_ErrMsg(string message) {
  string msg = message;
  msg += "\nProblem occurs at MC Counter " + to_string(Para.Counter) + ":\n";
  msg += DebugInfo(*Var.CurrGroup);
  return msg;
}

template <typename... TS>
std::string weight::ERR(std::string format, TS... args) {
  string msg = fmt::format(format, args...);
  msg += fmt::format("\nProblem occurs at MC Counter {}\n", Para.Counter);
  msg += DebugInfo(*Var.CurrGroup);
  return msg;
}

int weight::DynamicTest() {
  LOG_INFO("Start Dynamic Test...");
  LOG_INFO(DebugInfo(*Var.CurrGroup));
  //=================== Tau variable check ===================================//
  // for (int i = 1; i < Var.CurrGroup->TauNum / 2; i++)
  //   ASSERT_ALLWAYS(Equal(Var.Tau[2 * i], Var.Tau[2 * i + 1], 1.e-10),
  //                  ERR("Odd and Even are not the same! {0} vs {1}",
  //                      Var.Tau[2 * i], Var.Tau[2 * i + 1]));
  //================== External Variable
  //=========================================//
  // ASSERT_ALLWAYS(Equal(Var.Tau[0], 0.0, 1.0e-10), ERR("Tau 0 is not zero!"));
  ASSERT_ALLWAYS(Equal(Para.ExtMomTable[Var.CurrExtMomBin].data(),
                       Var.LoopMom[0].data(), D, 1.0e-8),
                 ERR("ExtMom is inconsistent! Bin: {0}; Mom: {1} vs {2}\n",
                     Var.CurrExtMomBin, Para.ExtMomTable[Var.CurrExtMomBin][0],
                     Var.LoopMom[0][0]));

  // ASSERT_ALLWAYS(
  //     Equal(Var.LoopMom[1].norm(), Para.Kf, 1.0e-8),
  //     ERR("ExtLegMom 1 is not on the Fermi surface! {0} not on Kf: {1}\n",
  //         Var.LoopMom[1].norm(), Para.Kf));

  // ASSERT_ALLWAYS(
  //     Equal(Var.LoopMom[2].norm(), Para.Kf, 1.0e-8),
  //     ERR("ExtLegMom 2 is not on the Fermi surface! {0} not on Kf: {1}\n",
  //         Var.LoopMom[2].norm(), Para.Kf));

  //=================== check NaN and Excited in weight ======================//
  for (auto &group : Groups) {
    // for (auto &diag : group.Diag) {
    //   for (auto i = 0; i < group.GNum; i++) {
    //     ASSERT_ALLWAYS(!std::isnan(diag.G[i]->Weight),
    //                    ERR("G weight is a NaN!\n Diag: {0} in Group: {1}\n",
    //                        diag.ID, group.ID));
    //     ASSERT_ALLWAYS(diag.G[i]->Excited == false,
    //                    ERR("One G is Excited!\n Diag: {0} in Group: {1}\n",
    //                        diag.ID, group.ID));
    //   }
    //   for (auto i = 0; i < group.Ver4Num; i++) {
    //     ASSERT_ALLWAYS(
    //         !std::isnan(diag.Ver4[i]->Weight[0] || diag.Ver4[i]->Weight[1]),
    //         ERR("One 4-Ver is NaN!\n Diag: {0} in Group: {1}\n", diag.ID,
    //             group.ID));
    //     ASSERT_ALLWAYS(diag.Ver4[i]->Excited[0] == false ||
    //                        diag.Ver4[i]->Excited[1] == false,
    //                    ERR("One 4-Ver is Excited!\n Diag: {0} in Group:
    //                    {1}\n",
    //                        diag.ID, group.ID));
    //   }
    // check diagram weight
    //   ASSERT_ALLWAYS(!std::isnan(diag.Weight),
    //                  ERR("Diag {} weight is a NaN!\n", diag.ID));
    // }
    // check group weight
    ASSERT_ALLWAYS(!std::isnan(group.Weight),
                   ERR("Group {} weight is a NaN!\n", group.ID));
  }

  //===============  Test if the weight is reproducible ============//
  // ChangeGroup(*Var.CurrGroup, true); // force to recalculate newweight
  double Weight = GetNewWeight(*Var.CurrGroup);

  // for (auto &diag : Var.CurrGroup->Diag) {
  //   for (auto i = 0; i < Var.CurrGroup->GNum; i++) {
  //     if (!Equal(diag.G[i]->NewWeight, diag.G[i]->Weight, 1.e-8)) {
  //       ASSERT_ALLWAYS(Equal(diag.G[i]->NewWeight, diag.G[i]->Weight, 1.e-8),
  //                      ERR("G Weight is different: {0} vs {1}\n",
  //                          diag.G[i]->NewWeight, diag.G[i]->Weight));
  //     }
  //     for (auto i = 0; i < Var.CurrGroup->Ver4Num; i++) {
  //       if (abs(diag.Ver4[i]->NewWeight[0] - diag.Ver4[i]->Weight[0])
  //       > 1.e-8)
  //         cout << "Error" << diag.Ver4[i]->NewWeight[0] << " vs "
  //              << diag.Ver4[i]->Weight[0] << endl;
  //       if (abs(diag.Ver4[i]->NewWeight[1] - diag.Ver4[i]->Weight[1])
  //       > 1.e-8)
  //         cout << "Error" << diag.Ver4[i]->NewWeight[1] << " vs "
  //              << diag.Ver4[i]->Weight[1] << endl;

  //       // ASSERT_ALLWAYS(
  //       //     Equal(diag.Ver4[i]->NewWeight[0],
  //       //     diag.Ver4[i]->Weight[0], 1.e-8), ERR("Ver4 Weight is
  //       different:
  //       //     {0} vs {1}\n",
  //       //         diag.Ver4[i]->NewWeight[0], diag.Ver4[i]->Weight[0]));
  //       // ASSERT_ALLWAYS(
  //       //     Equal(diag.Ver4[i]->NewWeight[1],
  //       //     diag.Ver4[i]->Weight[1], 1.e-8), ERR("Ver4 Weight is
  //       different:
  //       //     {0} vs {1}\n",
  //       //         diag.Ver4[i]->NewWeight[1], diag.Ver4[i]->Weight[1]));
  //     }
  //   }
  // }

  // if (abs(Weight - Var.CurrGroup->Weight) > 1.e-8)
  //   cout << "Error" << Weight << " vs " << Var.CurrGroup->Weight << endl;

  // ASSERT_ALLWAYS(
  //     Equal(Weight, Var.CurrGroup->Weight, 1.e-8),
  //     ERR("Weight is different: {0} vs {1}\n", Weight,
  //     Var.CurrGroup->Weight));

  // don't forget apply changes so that all Excited set to false
  // RejectChange(*Var.CurrGroup);

  //====== check reducibility ============================//
  // for (auto &diag : Var.CurrGroup->Diag) {
  //   for (int i = 0; i < Var.CurrGroup->Ver4Num; i++) {
  //     if (Para.UseVer4 == false) {
  //       vertex &Ver = *(diag.Ver[i]);
  //       if (IsInteractionReducible(Ver.LoopBasis[0], Var.CurrGroup->LoopNum))
  //         ASSERT_ALLWAYS(
  //             Equal(Ver.Weight[0], 0.0, 1.0e-10),
  //             ERR("Ver {0} is reducible but with finite Weight = {1}\n", i,
  //                 Ver.Weight[0]));

  //       if (IsInteractionReducible(Ver.LoopBasis[1], Var.CurrGroup->LoopNum))
  //         ASSERT_ALLWAYS(
  //             Equal(Ver.Weight[1], 0.0, 1.0e-10),
  //             ERR("Ver {0} is reducible but with finite Weight = {1}\n", i,
  //                 Ver.Weight[1]));
  //     }
  //   }
  // }

  LOG_INFO("Dynamic check past!");
  return 0;
}

template <typename T> bool IsInPool(const T *Pointer, const T *Pool, int Num) {
  if (Pointer < Pool && Pointer >= &Pool[Num])
    return false;
  else
    return true;
}

int weight::StaticTest() {

  // momentum InL = (sin(1.5), cos(1.5)) * Para.Kf;
  // momentum InR = (sin(2.5), cos(2.5)) * Para.Kf;
  // momentum DirTran = Para.ExtMomTable[ExtMomBinSize / 2];
  // Var.Tau[0] = Para.Beta / 3.0;
  // Var.Tau[1] = Para.Beta / 2.0;
  // int DiagIndex = 0;
  // int Level = 0;
  // int LoopNum = 1;
  // int LVerLoopNum = 0;
  // int TauIndex = 0;
  // int LoopIndex = 0;
  // double Tau=Var.Tau[1]-Var.Tau[0];

  // double GL2R = Fermi.Green(Tau, Internal2, UP, 0, Var.CurrScale);
  // double GR2L = Fermi.Green(TauR2L, Internal, UP, 0, Var.CurrScale);

  // int DiagIndex = OneLoop(InL, InR, DirTran, LoopNum, LVerLoopNum, TauIndex,
  //                         LoopIndex, DiagIndex, Level, T, false, false);
  // ASSERT_ALLWAYS(abs(_Weight[Level][DiagIndex][0]-)

  // check if all pointer is initilized
  // for (auto &group : Groups) {
  //   for (auto &dig : group.Diag) {
  //     for (auto i = 0; i < group.GNum; i++)
  //       ASSERT_ALLWAYS(
  //           dig.G[i] != nullptr,
  //           fmt::format("G {0} in Group {1} is null!\n", i, group.ID));
  //     for (auto i = 0; i < group.Ver4Num; i++) {
  //       ASSERT_ALLWAYS(
  //           dig.Ver4[i] != nullptr,
  //           fmt::format("Ver4 {0} in Group {1} is null!\n", i, group.ID));
  //     }
  //   }
  // }

  // check if pointers are corrected pointed to the object in Pool
  // ASSERT_ALLWAYS(Groups[0].Diag[0].G[0] == &Pool.GPool[0],
  //                "The first G in the first diagram does not pointed to the "
  //                "first G in GPool!");
  // if (Groups[0].Ver4Num > 0) {
  //   ASSERT_ALLWAYS(
  //       Groups[0].Diag[0].Ver4[0] == &Pool.Ver4Pool[0],
  //       "The first Ver4 in the first diagram does not pointed to the "
  //       "first Ver4 in Ver4Pool!");
  // }

  // for (auto &group : Groups) {
  //   for (auto &dig : group.Diag) {
  //     for (auto i = 0; i < group.GNum; i++)
  //       ASSERT_ALLWAYS(IsInPool(dig.G[i], Pool.GPool.data(), group.GNum),
  //                      fmt::format("G {0} in Group {1} is not in the
  //                      GPool\n!",
  //                                  i, group.ID));
  //     for (auto i = 0; i < group.Ver4Num; i++) {
  //       ASSERT_ALLWAYS(
  //           IsInPool(dig.Ver4[i], Pool.Ver4Pool.data(), group.Ver4Num),
  //           fmt::format("Ver4 {0} in Group {1} is not in the Ver4Pool\n!", i,
  //                       group.ID));
  //     }
  //   }
  // }

  return 0;
}
