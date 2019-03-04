#ifndef markov_H
#define markov_H

#include "global.h"
#include "utility/rng.h"
#include "weight.h"
#include <string>
#include <unordered_map>
#include <vector>

namespace mc {
using namespace std;
const int MCUpdates = 7;

typedef array<double, ExtMomBinSize> polar;

class markov {
public:
  markov();
  long long Counter;

  void PrintMCInfo();
  void PrintDeBugMCInfo();
  void AdjustGroupReWeight();

  // MC updates
  void ChangeTau();
  void ChangeMomentum();
  void IncreaseOrder();
  void DecreaseOrder();
  void ChangeGroup();

  // for RG only
  void RotateExtMom();
  void ChangeScale();

  void Measure();
  void SaveToFile();

  int DynamicTest();

  // MC variables
  diag::weight Weight;
  diag::variable &Var;
  vector<diag::group> &Groups;

private:
  // polarizatoin for each group
  unordered_map<int, polar> Polar;

  // polarizatoin for each group at the zero momentumr;
  unordered_map<int, double> PolarStatic;

  // beta function for 4-vertex
  double Ver4Flow[ScaleBinSize][InInAngBinSize][InOutAngBinSize];
  double Partition;

  // MC updates

  double ShiftExtK(const int &, int &);
  double ShiftK(const momentum &, momentum &, int Scale);
  double ShiftTau(const double &, double &);

  double GetNewExtK(momentum &, int Scale);
  double GetNewTau(double &);
  double GetNewK(momentum &, int Scale);
  double RemoveOldTau(double &);
  double RemoveOldK(momentum &, int Scale);

  // MC updates information
  std::string UpdatesName[MCUpdates];
  double Accepted[MCUpdates][MaxGroupNum];
  double Proposed[MCUpdates][MaxGroupNum];

  enum Updates {
    INCREASE_ORDER = 0,
    DECREASE_ORDER,
    CHANGE_GROUP,
    CHANGE_TAU,
    CHANGE_MOM,
    ROTATE_EXTK,
    CHANGE_SCALE,
    END
  };
  std::string _DetailBalanceStr(Updates op);
};
}; // namespace mc

#endif
