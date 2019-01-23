#ifndef markov_H
#define markov_H

#include "global.h"
#include "utility/rng.h"
namespace mc
{
const int MCUpdates = 3;
class Markov
{
public:
  long long *Counter;
  double Beta;
  int Order;
  // weight::GClass* G;
  // weight::WClass* W;
  RandomFactory *RNG;

  void Hop(int);
  void PrintDetailBalanceInfo();

  void ChangeTau();
  void ChangeMomentum();
  void ChangeGroup();

private:
  double ProbofCall[MCUpdates];
  double SumofProbofCall[MCUpdates];
  std::string OperationName[MCUpdates];
  double Accepted[MCUpdates][MaxOrder];
  double Proposed[MCUpdates][MaxOrder];

  double RandomPickK();
  double RandomPickTau();
  enum Operations
  {
    CHANGE_TAU = 0,
    CHANGE_MOMENTUM,
    CHNAGE_GROUP,
    END
  };
  std::string _DetailBalanceStr(Operations op);
  std::string _CheckBalance(Operations op1, Operations op2);
};
}; // namespace mc

#endif
