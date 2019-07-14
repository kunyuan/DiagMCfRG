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

  msg += string(80, '=') + "\n";
  msg += "LoopMom: \n";
  for (int i = 0; i < Group.LoopNum; i++)
    msg +=
        ToString(Var.LoopMom[i]) + "=" + ToString(Var.LoopMom[i].norm()) + "\n";

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

  for (int scale = 0; scale < ScaleBinSize; ++scale)
    for (int angle = 0; angle < AngBinSize; ++angle)
      for (int q = 0; q < ExtMomBinSize; ++q) {
        ASSERT_ALLWAYS(!std::isnan(VerQTheta.EffInteraction[scale][angle][q]),
                       ERR("VerQTheta contains a NaN! Index: {0}, {1}, {2}\n",
                           scale, angle, q));
      }

  //=================== check NaN and Excited in weight ======================//
  for (auto &group : Groups) {
    ASSERT_ALLWAYS(!std::isnan(real(group.Weight)),
                   ERR("Group {} weight is a NaN!\n", group.ID));
  }

  for (auto &group : Groups) {
    ASSERT_ALLWAYS(!std::isnan(imag(group.Weight)),
                   ERR("Group {} weight is a NaN!\n", group.ID));
  }

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
