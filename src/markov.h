#ifndef markov_H
#define markov_H

#include "global.h"
#include "utility/estimator.h"
#include "utility/rng.h"
#include "weight.h"
#include <string>
#include <vector>

namespace mc {
using namespace std;
const int MCUpdates = 3;

typedef array<double, ExtMomBinSize> polar;

class markov {
public:
  long long Counter;

  void Initialization(std::string FilePrefix);
  void Hop(const int);
  void PrintMCInfo();
  void AdjustGroupReWeight();

  void Measure();
  void SaveToFile(std::string FilePrefix);

  int DynamicTest();

private:
  diag::weight Weight;
  // double SumofProbofCall[MCUpdates];
  std::string UpdateName[MCUpdates];
  double Accepted[MCUpdates][MaxGroupNum];
  double Proposed[MCUpdates][MaxGroupNum];

  void ChangeTau();
  void ChangeMomentum();
  void ChangeGroup();

  int RandomPickExtK(const int &, double &);
  void RandomPickK(const momentum &, momentum &, double &);
  double RandomPickTau(const double &, double &);
  enum Updates { CHANGE_GROUP = 0, CHANGE_TAU, CHANGE_MOM, END };
  std::string _DetailBalanceStr(Updates op);

  // polarizatoin for each group
  vector<polar> Polar;
  // polarizatoin for each group at the zero momentumr;
  EstimatorBundle<double> PolarStatic;
};
}; // namespace mc

#endif
