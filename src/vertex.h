#ifndef vertex_H
#define vertex_H

// #include "utility/vector.h"
#include "global.h"

double sum2(const momentum &);
double norm2(const momentum &);

namespace diag {
double Green(double Tau, const momentum &Momentum, spin Spin, int GType);
double Interaction(double Tau, const momentum &Momentum, int VerType);
} // namespace diag
#endif
