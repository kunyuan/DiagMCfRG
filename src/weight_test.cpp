#include "weight.h"

using namespace diag;
using namespace vertex;
using namespace std;

string weight::_DebugInfo() {
  string msg;
  msg = string(80, '=') + "\n";
  msg += "\nMC Counter: " + to_string(Para.Counter) + ":\n";
  msg += "Current Group Info:\n " + ToString(*Var.CurrGroup);

  msg += string(80, '=') + "\n";
  msg += "GWeight: \n";
  for (int i = 0; i < Var.CurrGroup->GNum; i++)
    msg += ToString(Var.CurrGroup->Diag[0].G[i]->Weight) + "; ";
  msg += "\n";

  msg += "NewGWeight: \n";
  for (int i = 0; i < Var.CurrGroup->GNum; i++)
    msg += ToString(Var.CurrGroup->Diag[0].G[i]->NewWeight) + "; ";
  msg += "\n";

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

  msg += "VerWeight: \n";
  for (int i = 0; i < Var.CurrGroup->Ver4Num; i++) {
    msg += ToString(Var.CurrGroup->Diag[0].Ver[i]->Weight[0]) + ", ";
    msg += ToString(Var.CurrGroup->Diag[0].Ver[i]->Weight[1]) + "; ";
  }
  msg += "\n";

  msg += string(80, '=') + "\n";
  msg += "LoopMom: \n";
  for (int d = 0; d < D; d++) {
    for (int i = 0; i < Var.CurrGroup->LoopNum; i++)
      msg += ToString(Var.LoopMom[i][d]) + ", ";
    msg += "\n";
  }
  msg += "\n";

  msg += string(80, '=') + "\n";
  msg += "Tau: \n";
  for (int i = 0; i < Var.CurrGroup->TauNum; i++)
    msg += ToString(Var.Tau[i]) + ", ";

  msg += "\n";
  return msg;
}

string weight::_ErrMsg(string message) {
  string msg = message;
  msg += "\nProblem occurs at MC Counter " + to_string(Para.Counter) + ":\n";
  msg += _DebugInfo();
  return msg;
}

int weight::DynamicTest() {
  LOG_INFO("Start Dynamic Test...");
  LOG_INFO(_DebugInfo());
  //=================== Tau variable check ===================================//
  for (int i = 1; i < Var.CurrGroup->TauNum / 2; i++)
    ASSERT_ALLWAYS(Equal(Var.Tau[2 * i], Var.Tau[2 * i + 1], 1.e-10),
                   "Odd and Even are not the same! " << Var.Tau[2 * i] << " vs "
                                                     << Var.Tau[2 * i + 1]);
  //================== External Variable
  //=========================================//
  ASSERT_ALLWAYS(Equal(Var.Tau[0], 0.0, 1.0e-10),
                 _ErrMsg("Tau 0 is not zero!"));
  ASSERT_ALLWAYS(
      Equal(Var.ExtMomTable[Var.CurrExtMomBin].data(), Var.LoopMom[0].data(), D,
            1.0e-8),
      _ErrMsg("ExtMom is inconsistent! Bin: " + ToString(Var.CurrExtMomBin) +
              "\n" + ToString(Var.ExtMomTable[Var.CurrExtMomBin][0]) + " vs " +
              ToString(Var.LoopMom[0][0])));

  //=================== check NaN and Excited in weight ======================//
  for (auto &group : Groups) {
    for (auto &diag : group.Diag) {
      for (auto i = 0; i < group.GNum; i++) {
        ASSERT_ALLWAYS(
            !isnan(diag.G[i]->Weight),
            _ErrMsg("One G weight is a NaN!\n Diag: " + ToString(diag) +
                    "\n Group: " + ToString(group)));
        ASSERT_ALLWAYS(diag.G[i]->Excited == false,
                       _ErrMsg("One G is Excited!\n Diag: " + ToString(diag) +
                               "\n Group: " + ToString(group)));
      }
      for (auto i = 0; i < group.Ver4Num; i++) {
        if (UseVertex4) {
          ASSERT_ALLWAYS(
              !isnan(diag.Ver4[i]->Weight),
              _ErrMsg("One 4-Ver weight is a NaN!\n Diag: " + ToString(diag) +
                      "\n Group: " + ToString(group)));
          ASSERT_ALLWAYS(
              diag.Ver4[i]->Excited == false,
              _ErrMsg("One Ver4 is Excited!\n Diag: " + ToString(diag) +
                      "\n Group: " + ToString(group)));
        } else {
          ASSERT_ALLWAYS(
              !isnan(diag.Ver[i]->Weight[0] || diag.Ver[i]->Weight[1]),
              _ErrMsg("One Ver weight is a NaN!\n Diag: " + ToString(diag) +
                      "\n Group: " + ToString(group)));
          ASSERT_ALLWAYS(
              diag.Ver[i]->Excited[0] == false ||
                  diag.Ver[i]->Excited[1] == false,
              _ErrMsg("One Ver is Excited!\n Diag: " + ToString(diag) +
                      "\n Group: " + ToString(group)));
        }
      }
      // check diagram weight
      ASSERT_ALLWAYS(
          !isnan(diag.Weight),
          _ErrMsg("Diag weight is a NaN!\n Diag: " + ToString(diag)));
    }
    // check group weight
    ASSERT_ALLWAYS(
        !isnan(group.Weight),
        _ErrMsg("Group weight is a NaN!\n Diag: " + ToString(group)));
  }

  //===============  Test if the weight is reproducible ============//
  ChangeGroup(*Var.CurrGroup, true); // force to recalculate newweight
  double Weight = GetNewWeight(*Var.CurrGroup);

  for (auto &diag : Var.CurrGroup->Diag) {
    for (auto i = 0; i < Var.CurrGroup->GNum; i++) {
      if (!Equal(diag.G[i]->NewWeight, diag.G[i]->Weight, 1.e-8)) {
        ASSERT_ALLWAYS(
            Equal(diag.G[i]->NewWeight, diag.G[i]->Weight, 1.e-8),
            _ErrMsg("G Weight is different: " + ToString(diag.G[i]->NewWeight) +
                    " vs " + ToString(diag.G[i]->Weight)));
      }
      for (auto i = 0; i < Var.CurrGroup->Ver4Num; i++) {
        if (UseVertex4) {
          ASSERT_ALLWAYS(
              Equal(diag.Ver4[i]->NewWeight, diag.Ver4[i]->Weight, 1.e-8),
              _ErrMsg("Ver4 Weight is different: " +
                      ToString(diag.Ver4[i]->NewWeight) + " vs " +
                      ToString(diag.Ver4[i]->Weight)));
        } else {
          ASSERT_ALLWAYS(
              Equal(diag.Ver[i]->NewWeight[0], diag.Ver[i]->Weight[0], 1.e-8),
              _ErrMsg("Ver Weight is different: " +
                      ToString(diag.Ver[i]->NewWeight[0]) + " vs " +
                      ToString(diag.Ver[i]->Weight[0])));
          ASSERT_ALLWAYS(
              Equal(diag.Ver[i]->NewWeight[1], diag.Ver[i]->Weight[1], 1.e-8),
              _ErrMsg("Ver Weight is different: " +
                      ToString(diag.Ver[i]->NewWeight[1]) + " vs " +
                      ToString(diag.Ver[i]->Weight[1])));
        }
      }
    }
  }

  ASSERT_ALLWAYS(Equal(Weight, Var.CurrGroup->Weight, 1.e-8),
                 _ErrMsg("Weight is different: " + ToString(Weight) + " vs " +
                         ToString(Var.CurrGroup->Weight)));
  // don't forget apply changes so that all Excited set to false
  RejectChange(*Var.CurrGroup);

  //====== check reducibility ============================//
  for (auto &diag : Var.CurrGroup->Diag) {
    for (int i = 0; i < Var.CurrGroup->Ver4Num; i++) {
      vertex &Ver = *(diag.Ver[i]);
      if (IsInteractionReducible(Ver.LoopBasis[0], Var.CurrGroup->LoopNum)) {
        ASSERT_ALLWAYS(Equal(Ver.Weight[0], 0.0, 1.0e-10),
                       _ErrMsg("Ver " + ToString(i) +
                               " is reducible but with finite Weight=" +
                               ToString(Ver.Weight[0])));
      }
      if (IsInteractionReducible(Ver.LoopBasis[1], Var.CurrGroup->LoopNum)) {
        ASSERT_ALLWAYS(Equal(Ver.Weight[1], 0.0, 1.0e-10),
                       _ErrMsg("Ver " + ToString(i) +
                               " is reducible but with finite Weight=" +
                               ToString(Ver.Weight[1])));
      }
    }
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
  for (auto &group : Groups) {
    for (auto &dig : group.Diag) {
      for (auto i = 0; i < group.GNum; i++)
        ASSERT_ALLWAYS(dig.G[i] != nullptr,
                       "G " << i << " in Group " << group.ID << " is null!");
      for (auto i = 0; i < group.Ver4Num; i++) {
        ASSERT_ALLWAYS(dig.Ver[i] != nullptr,
                       "Ver " << i << " in Group " << group.ID << " is null!");
        ASSERT_ALLWAYS(dig.Ver4[i] != nullptr,
                       "Ver4 " << i << " in Group " << group.ID << " is null!");
      }
    }
  }

  // check if pointers are corrected pointed to the object in Pool
  ASSERT_ALLWAYS(Groups[0].Diag[0].G[0] == &Pool.GPool[0],
                 "The first G in the first diagram does not pointed to the "
                 "first G in GPool!");
  if (Groups[0].Ver4Num > 0) {
    ASSERT_ALLWAYS(Groups[0].Diag[0].Ver[0] == &Pool.VerPool[0],
                   "The first Ver in the first diagram does not pointed to the "
                   "first Ver in VerPool!");
    ASSERT_ALLWAYS(
        Groups[0].Diag[0].Ver4[0] == &Pool.Ver4Pool[0],
        "The first Ver4 in the first diagram does not pointed to the "
        "first Ver4 in Ver4Pool!");
  }

  for (auto &group : Groups) {
    for (auto &dig : group.Diag) {
      for (auto i = 0; i < group.GNum; i++)
        ASSERT_ALLWAYS(IsInPool(dig.G[i], Pool.GPool.data(), group.GNum),
                       "G " << i << " in Group " << group.ID
                            << " is not in the GPool!");
      for (auto i = 0; i < group.Ver4Num; i++) {
        ASSERT_ALLWAYS(IsInPool(dig.Ver[i], Pool.VerPool.data(), group.Ver4Num),
                       "Ver " << i << " in Group " << group.ID
                              << " is not in the VerPool!");
        ASSERT_ALLWAYS(
            IsInPool(dig.Ver4[i], Pool.Ver4Pool.data(), group.Ver4Num),
            "Ver4 " << i << " in Group " << group.ID
                    << " is not in the Ver4Pool!");
      }
    }
  }

  return 0;
}
