#ifndef vertex_H
#define vertex_H

// #include "utility/vector.h"
#include "global.h"
#include <array>

double sum2(const momentum &);
double norm2(const momentum &);

namespace diag {

const int MAXSIGMABIN = 100000;

class fermi {
public:
  fermi();
  double Green(double Tau, const momentum &Momentum, spin Spin, int GType,
               double Scale = 0);

private:
  // beyond which the expense sigma function will be called
  double UpperBound, LowerBound;
  double DeltaK;
  double UpperBound2, LowerBound2; // lower upbound for better sigma
  double DeltaK2;
  double PhyGreen(double Tau, const momentum &Mom, int GType, double Scale = 0);
  double FockSigma(const momentum &Mom);
  double BuildFockSigma();
  double Fock(double k);
  // warning: this function only works for T=0 and 3D!!!!
  double GetSigma(double k);
  double Sigma[MAXSIGMABIN];
  double Sigma2[MAXSIGMABIN];
};

class verQTheta {
public:
  verQTheta();
  double Interaction(const momentum &InL, const momentum &InR,
                     const momentum &Transfer, double Tau, int VerType);

  void Measure(const momentum &InL, const momentum &InR, const int QIndex,
               int Order, double *Tau, double *Weight, double WeightFactor);
  void Update(double Ratio = 1.0);
  void Save();
  void ClearStatis();
  void ResetIRScale(int IRScaleBin);
  double *EffInteraction;
  double *DiffInteraction;
  double *IntInteraction;

  double &EffInter(int Angle, int ExtQ, int Tau);
  double &DiffInter(int Order, int Angle, int ExtQ, int Tau);
  // double EffInteraction[ScaleBinSize + 1][AngBinSize][ExtMomBinSize];
  // double DiffInteraction[MaxOrder][ScaleBinSize +
  // 1][AngBinSize][ExtMomBinSize]; double IntInteraction[MaxOrder][ScaleBinSize
  // + 1][AngBinSize][ExtMomBinSize];

  // double TauBasis[TauBinSize][TauBasisNum];

  double Normalization;
  double PhyWeight;

  int QIndex;
  int AngleIndex;
  int OrderIndex;
};

class verQ {
public:
  double Interaction(double Tau, const momentum &Momentum, int VerType,
                     double Scale = 0);
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

double Index2Scale(const int &Index);
int Scale2Index(const double &Scale);

double Index2Tau(const int &Index);
int Tau2Index(const double &Tau);

}; // namespace diag
#endif
