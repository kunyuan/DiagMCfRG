#ifndef vertex_H
#define vertex_H

// #include "utility/vector.h"
#include "global.h"

double sum2(const momentum &);
double norm2(const momentum &);

namespace diag {

const int MAXSIGMABIN = 8192;
class fermi {
public:
  double Green(double Tau, const momentum &Momentum, spin Spin, int GType);

private:
  // beyond which the expense sigma function will be called
  const double UpBound = 10.0;
  double PhyGreen(double Tau, const momentum &Mom);
  double FockSigma(const momentum &Mom);
  double BuildFockSigma();
  double Fock(double k);
  double GetSigma(double k);
  double Sigma[MAXSIGMABIN];
};

class bose {
public:
  double Interaction(double Tau, const momentum &Momentum, int VerType);
};

}; // namespace diag
#endif
