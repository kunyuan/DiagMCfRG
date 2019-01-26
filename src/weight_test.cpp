#include "weight.h"

using namespace diag;
using namespace vertex;
using namespace std;

string weight::_ErrMsg(string message) {
  string msg = message;
  msg += "\n Problem occurs at Step " + to_string(Para.Counter) + ":\n";
  msg += "Current Group Info:\n" + ToString(*Var.CurrGroup);
  msg += "Current Weight: " + to_string(Var.CurrWeight) + "\n";
  return msg;
}

int weight::DynamicTest() {
  //=================== check NaN in weight  ======================//
  for (auto &group : Groups) {
    for (auto &diag : group.Diag) {
      for (auto &g : diag.G) {
        ASSERT_ALLWAYS(
            !isnan(g->Weight),
            _ErrMsg("One G weight is a NaN!\n Diag: " + ToString(diag) +
                    "\n Group: " + ToString(group)));
      }
      if (UseVertex4) {
        for (auto &v : diag.Ver4) {
          ASSERT_ALLWAYS(
              !isnan(v->Weight),
              _ErrMsg("One 4-Ver weight is a NaN!\n Diag: " + ToString(diag) +
                      "\n Group: " + ToString(group)));
        }
      } else {
        for (auto &v : diag.Ver) {
          ASSERT_ALLWAYS(
              !isnan(v->Weight[0] || v->Weight[1]),
              _ErrMsg("One Ver weight is a NaN!\n Diag: " + ToString(diag) +
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
  // check total weight
  ASSERT_ALLWAYS(!isnan(Var.CurrWeight), _ErrMsg("Total weight is a NaN!\n"));

  //===============  Test if the weight is reproducible ============//
  ChangeGroup(*Var.CurrGroup, true); // force to recalculate newweight
  double Weight = GetNewWeight(*Var.CurrGroup);
  ASSERT_ALLWAYS(Equal(Weight, Var.CurrWeight, 1.e-6),
                 _ErrMsg("Weight is different: " + to_string(Weight) + " vs " +
                         to_string(Var.CurrWeight)));
  return 0;
}
