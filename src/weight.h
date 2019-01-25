#ifndef weight_H
#define weight_H

#include "diagram.h"
#include "vertex.h"
#include "utility/rng.h"
extern parameter Para;
extern RandomFactory Random;

namespace diag
{
using namespace std;
struct variable
{
    group *CurrGroup;
    long int CurrVersion;
    double CurrWeight;                   //current weight of the system
    array<momentum, MaxLoopNum> LoopMom; //all momentum loop variables
    array<double, MaxTauNum> Tau;        //all tau variables
    array<int, MaxLoopNum> LoopSpin;     //all spin variables
};

class weight
{
  public:
    // double ChangeTemperature(double NewBeta);
    double GetNewWeight(group &); //return the current weight
    void ChangeMom(group &, int Index);
    void ChangeTau(group &, int TauIndex);          //two tau index on the two sides of interaction
    void ChangeGroup(group &, bool Forced = false); //recalculate the weights in one group
    void AcceptChange(group &);
    void ReadDiagrams(std::string FilePrefix);
    void Initialization();

  private:
    variable Variable; //The variable of the integral
    vector<group> GroupList;
    pool Pool; //Pool to store indepdent G, Vertex, and 4-Vertex
    struct
    {
        int Num;
        array<double, MaxGNum> Weight;
        array<green *, MaxGNum> Index;
    } NewG;
    struct
    {
        int Num;
        array<double, MaxVer4Num * 2> Weight;
        array<vertex *, MaxVer4Num * 2> Index;
    } NewVer;
    struct
    {
        int Num;
        array<double, MaxVer4Num> Weight;
        array<vertex4 *, MaxVer4Num> Index;
    } NewVer4;
    momentum _Mom(const loop &LoopBasis, const int &LoopNum);
    array<array<double, diag::MaxBranchNum>, MaxVer4Num> _SpinCache;
};

}; // namespace diag

#endif
