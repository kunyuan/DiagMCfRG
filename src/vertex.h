#ifndef vertex_H
#define vertex_H

// #include "utility/vector.h"
#include "global.h"

double sum2(const momentum &);
double norm2(const momentum &);

namespace diag {

const int MAXSIGMABIN = 100000;

class fermi {
public:
  fermi();
  double Green(double Tau, const momentum &Momentum, spin Spin, int GType,
               int Scale = 0);

private:
  // beyond which the expense sigma function will be called
  double UpperBound, LowerBound;
  double DeltaK;
  double UpperBound2, LowerBound2; // lower upbound for better sigma
  double DeltaK2;
  double PhyGreen(double Tau, const momentum &Mom, int Scale = 0);
  double FockSigma(const momentum &Mom);
  double BuildFockSigma();
  double Fock(double k);
  // warning: this function only works for T=0 and 3D!!!!
  double GetSigma(double k);
  double Sigma[MAXSIGMABIN];
  double Sigma2[MAXSIGMABIN];
};

class bose {
public:
  bose();
  double Interaction(double Tau, const momentum &Momentum, int VerType,
                     int Scale = 0);

  double Interaction(const momentum &InL, const momentum &InR,
                     const momentum &Transfer, int VerType, int Scale = 0);

  double EffInteraction[ScaleBinSize][AngBinSize][ExtMomBinSize];
  double DiffInteraction[ScaleBinSize][AngBinSize][ExtMomBinSize];
  double Normalization[ScaleBinSize];
};

class verfunc {
public:
  verfunc();
  void Vertex4(const momentum &InL, const momentum &InR, const momentum &OutL,
               const momentum &OutR, int Ver4TypeDirect, int Ver4TypeExchange,
               double &Direct, double &Exchange);

  // private:
  //   double Ver4AtUV[InInAngBinSize][InOutAngBinSize];
  //   double Angle2D(const momentum &K1, const momentum &K2);
  //   double Index2Angle(const int &Index, const int &AngleNum);
  //   int Angle2Index(const double &Angle, const int &AngleNum);
  //   void _TestAngleIndex();
  //   void _TestAngle2D();
};

double Angle2D(const momentum &K1, const momentum &K2);
double Index2Angle(const int &Index, const int &AngleNum);
int Angle2Index(const double &Angle, const int &AngleNum);
void _TestAngleIndex();
void _TestAngle2D();

double Index2Mom(const int &Index);
int Mom2Index(const double &K);

}; // namespace diag
#endif
