#ifndef markov_H
#define markov_H

#include "utility/convention.h"
#include "utility/rng.h"
namespace mc {
const int MCUpdates = 3;
class Markov {
public:
    long long* Counter;
    double Beta;
    int Order;
    // weight::GClass* G;
    // weight::WClass* W;
    RandomFactory* RNG;

    void Hop(int);
    void PrintDetailBalanceInfo();

    void ChangeTau();
    void ChangeMomentum();
    void ChangeGroup();

private:
    double ProbofCall[MCUpdates];
    double SumofProbofCall[MCUpdates];
    std::string OperationName[MCUpdates];
    double Accepted[MCUpdates][MAX_ORDER];
    double Proposed[MCUpdates][MAX_ORDER];

    real RandomPickK();
    real RandomPickTau();
    enum Operations {
        CHANGE_TAU = 0,
        CHANGE_MOMENTUM,
        CHNAGE_GROUP,
        END
    };
    std::string _DetailBalanceStr(Operations op);
    std::string _CheckBalance(Operations op1, Operations op2);
}

#endif
