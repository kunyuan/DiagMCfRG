#ifndef weight_H
#define weight_H

#include "diagram.h"
#include "utility/rng.h"
#include "vertex.h"
extern parameter Para;
extern RandomFactory Random;

namespace diag {
using namespace std;
struct variable {
  group *CurrGroup;
  long int CurrVersion;
  int CurrExtMomBin; // current bin of the external momentum
  double CurrTau;    // current external tau
  double CurrFreq;   // current external tau
  double CurrScale;  // Current (Reference) Scale: Index=1, ..., ScaleBinSize
  int CurrIRScaleBin;
  array<momentum, MaxLoopNum> LoopMom; // all momentum loop variables
  // array<int, MaxLoopNum> LoopFreq;     // all momentum loop variables
  array<double, MaxTauNum> Tau;    // all tau variables
  array<int, MaxLoopNum> LoopSpin; // all spin variables
};

class weight {
public:
  vector<group> Groups;
  variable Var; // The variable of the integral

  // initialization, read diagrams, then initialize variables
  void ReadDiagrams();

  // MC updates related operations
  // double ChangeTemperature(double NewBeta);
  void ChangeMom(group &, int Index);
  void ChangeTau(group &, int TauIndex);
  // two tau index on the two sides of interaction
  void ChangeGroup(group &, bool Forced = false);
  // recalculate the weights in one group
  cmplx GetNewWeight(group &); // return the current weight
  void AcceptChange(group &);
  void RejectChange(group &);

  cmplx Evaluate(group &);

  void Measure(double WeightFactor);
  void Update(double Ratio);
  void ClearStatis();
  void ResetIRScale(int IRScaleBin);
  void Save();

  // run test in MC updates
  int DynamicTest();

  // Test before MC
  int StaticTest();

  string DebugInfo(group &);

private:
  pool Pool; // Pool to store indepdent G, Vertex, and 4-Vertex
  struct {
    int Num;
    array<double, MaxGNum> Weight;
    array<green *, MaxGNum> Index;
  } NewG;
  struct {
    int Num;
    array<double, MaxVer4Num> Weight;
    array<vertex4 *, MaxVer4Num> Index;
  } NewVer4;
  array<array<double, MaxBranchNum>, MaxVer4Num> _SpinCache;
  string _ErrMsg(string);

  void Initialization();

  void GetMom(const loop &LoopBasis, const int &LoopNum, momentum &Mom);
  momentum _Mom;
  momentum _InL;
  momentum _InR;
  momentum _OutL;
  momentum _OutR;

  // the spin cache to calculate vertex weight
  double _Tree[MaxOrder][MaxBranchNum];
  bool IsInteractionReducible(loop &, int LoopNum);
  bool IsInteractionReducible(loop &, loop &, int LoopNum);

  template <typename... TS> string ERR(string format, TS... args);

  fermi Fermi;
  verQ VerQ;
  verQTheta VerQTheta;
  verfunc VerFunc;

  cmplx fRG(int LoopNum);
  cmplx Ver4Loop(const momentum &InL, const momentum &InR,
                 const momentum &DirTran, const momentum &ExTran, int Order,
                 int LoopIndex, int Channel, int Type, int &Level,
                 int &DiagNum);
  cmplx Ver4Loop0(const momentum &InL, const momentum &InR,
                  const momentum &DirTran, const momentum &ExTran,
                  int LoopIndex, int &Level, int &DiagNum);
  cmplx Ver4Loop1(const momentum &InL, const momentum &InR,
                  const momentum &DirTran, const momentum &ExTran, int Order,
                  int LoopIndex, int Channel, int Type, int &Level,
                  int &DiagNum);

  cmplx _Weight[MaxOrder][MaxDiagNum][2];
  double _ExtTau[MaxOrder][MaxDiagNum][4];
  // double Ver4Loop2();
  // double Ver6Loop1();
};

}; // namespace diag

#endif
