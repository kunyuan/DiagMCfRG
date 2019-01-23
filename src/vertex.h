#ifndef vertex_H
#define vertex_H

#include "utility/vector.h"
#include "global.h"

namespace vertex
{
double Green(double Tau, const mom &Momentum, spin Spin, int GType);
double Interaction(double Tau, const mom &Momentum, spin Spin, int VerType);
} // namespace vertex
#endif
