#ifndef weight_H
#define weight_H

#include "diagram.h"
#include "utility/rng.h"
#include "vertex.h"
#include <vector>
extern parameter Para;
extern RandomFactory Random;

namespace diag {
using namespace std;
struct variable {
  group *CurrGroup;
  long int CurrVersion;
  int CurrExtMomBin; // current bin of the external momentum
  double CurrTau;    // current external tau
  double CurrScale;  // Current (Reference) Scale: Index=1, ..., ScaleBinSize
  int CurrIRScaleBin;
  double CurrWeight[MaxTauNum];
  array<momentum, MaxLoopNum> LoopMom; // all momentum loop variables
  array<double, MaxTauNum> Tau;        // all tau variables
  array<int, MaxLoopNum> LoopSpin;     // all spin variables
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
  double GetNewWeight(group &); // return the current weight
  void AcceptChange(group &);
  void RejectChange(group &);

  double Evaluate(group &);

  void Measure(double WeightFactor);
  void Update(double Ratio);
  void ClearStatis();
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

  double fRG(int LoopNum, int ID);
  int Vertex4(
      const momentum &InL, const momentum &InR, const momentum &DirTran,
      int LoopNum, int TauIndex, int LoopIndex, int Level,
      bool *Channel,     // three flags, calculate t, u, s or not
      int VerType = -1,  // -1: normal, 0: left(to project), 1: right(to diff)
      int LVerOrder = -1 // order of left vertex
  );

  int Bubble(
      const momentum &InL, const momentum &InR, const momentum &DirTran,
      int LoopNum, int TauIndex, int LoopIndex, int Level,
      bool *Channel,      // three flags, calculate t, u, s or not
      int VerType = -1,   // -1: normal, 0: left(to project), 1: right(to diff)
      int LVerOrder = -1, // order of left vertex
      bool IsPenguin = false);

  int Ver4Loop0(const momentum &InL, const momentum &InR,
                const momentum &DirTran, int TauIndex, int LoopIndex, int Level,
                int Type = 0 // 0: renormalized interaction, -2: bare coupling
  );

  int OneLoop(const momentum &InL, const momentum &InR, const momentum &DirTran,
              int LoopNum, int LVerLoopNum, int TauIndex, int LoopIndex,
              int Level,
              bool *Channel, // three flags, calculate t, u, s or not
              bool IsProjected = false, bool IsPenguin = false);

  int _DiagIndex[MaxLevel];
  int _ExtTau[MaxOrder][MaxLevel][MaxDiagNum][4];
  double _Weight[MaxLevel][MaxDiagNum][2];
  int _DiagNum;
  int _GlobalOrder;
  double _GL2R[MaxTauNum][MaxTauNum];
  double _GR2L[MaxTauNum][MaxTauNum];
  // double Ver4Loop2();
  // double Ver6Loop1();
  bool ALL[3] = {true, true, true};
  bool US[3] = {false, true, true};
  bool UT[3] = {true, true, false};
  bool ST[3] = {true, false, true};
  bool T[3] = {true, false, false};
  bool U[3] = {false, true, false};
  bool S[3] = {false, false, true};
  int LEFT = 0;
  int RIGHT = 1;
};

}; // namespace diag

#endif
