#include "vertex.h"
#include "global.h"
#include "utility/abort.h"
#include <cmath>
#include <iostream>

using namespace vertex;
using namespace std;

extern parameter Para;

double sum2(const momentum &Mom) {
  double Sum2 = 0.0;
  for (int i = 0; i < D; i++)
    Sum2 += Mom[i] * Mom[i];
  return Sum2;
}

double norm2(const momentum &Mom) { return sqrt(sum2(Mom)); }

double vertex::Interaction(double Tau, const momentum &Mom, int VerType) {
  if (VerType < 0)
    ABORT("VerType can not be " << VerType);
  double interaction = 8.0 * PI / (sum2(Mom) + Para.Mass2);

  if (VerType > 0) {
    // the interaction contains counter-terms
    interaction *= pow(Para.Mass2 / (sum2(Mom) + Para.Mass2), VerType);
    interaction *= pow(-1, VerType);
  }
  return interaction;
}

double PhyGreen(double Tau, const momentum &Mom) {
  // if tau is exactly zero, set tau=0^-
  double green, Ek;
  if (Tau == 0.0) {
    return EPS;
  }

  double s = 1.0;
  if (Tau < 0.0) {
    Tau += Para.Beta;
    s = -s;
  } else if (Tau >= Para.Beta) {
    Tau -= Para.Beta;
    s = -s;
  }
  Ek = sum2(Mom); // bare propagator

  // Ek=SelfEnergy(Mom);   //Fock diagram dressed propagator

  //// enforce an UV cutoff for the Green's function ////////
  // if(Ek>8.0*EF) then
  //   PhyGreen=0.0
  //   return
  // endif

  double x = Para.Beta * (Ek - Para.Mu) / 2.0;
  double y = 2.0 * Tau / Para.Beta - 1.0;
  if (x > 100.0)
    green = exp(-x * (y + 1.0));
  else if (x < -100.0)
    green = exp(x * (1.0 - y));
  else
    green = exp(-x * y) / (2.0 * cosh(x));

  green *= s;

  // cout << "x: " << x << ", y: " << y << ", G: " << green << endl;
  // cout << "G: " << green << endl;

  if (std::isnan(green))
    ABORT("Green is too large!" << Tau << " " << Ek << " " << green);
  return green;
}

double vertex::Green(double Tau, const momentum &Mom, spin Spin, int GType) {
  double green;
  if (GType == 0) {
    green = PhyGreen(Tau, Mom);
    cout << "Tau=" << Tau << endl;
    cout << "Mom=" << Mom[0] * Mom[0] + Mom[1] * Mom[1] << endl;
    cout << green << endl;
  } else
    ABORT("GType " << GType << " has not yet been implemented!");
  // return FakeGreen(Tau, Mom);
  return green;
}