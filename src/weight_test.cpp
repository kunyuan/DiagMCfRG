#include "weight.h"

using namespace diag;
using namespace vertex;
using namespace std;

string weight::_ErrMsg(string message) {
  string msg = message;
  msg += "\nProblem occurs at MC Counter " + to_string(Para.Counter) + ":\n";
  msg += "Current Group Info:\n " + ToString(*Var.CurrGroup) + "\n";
  return msg;
}

int weight::DynamicTest() {
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
  ASSERT_ALLWAYS(Equal(Weight, Var.CurrGroup->Weight, 1.e-8),
                 _ErrMsg("Weight is different: " + ToString(Weight) + " vs " +
                         ToString(Var.CurrGroup->Weight)));

  // don't forget apply changes so that all Excited set to false
  AcceptChange(*Var.CurrGroup);
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
