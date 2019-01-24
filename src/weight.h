#ifndef weight_H
#define weight_H

#include "diagram.h"
extern parameter Para;

namespace diag
{
class weight
{
  public:
    // double ChangeTemperature(double NewBeta);
    double ChangeMom(int Index);
    double ChangeTau(int IndexIn, int IndexOut); //IndexIn may be equal to IndexOut
    double GetWeight();                          //return the current weight
    void ReadDiagrams(std::string FilePrefix);

  private:
    double CurrWeight; //current weight of the system
    struct
    {
        double LoopMom[MaxOrder + 1][D];
        double Tau[2 * MaxOrder];
        int LoopSpin[MaxOrder + 1];
    } variable;
    vector<group> GroupList;
    vector<green> GPool;
    vector<vertex> VerPool;
    vector<vertex4> Ver4Pool;
    // void _ReadOneGroup(ifstream &, group &);
    void _AddNewGToPool(diagram &, vector<int> &, vector<loop> &, vector<int> &);
    void _AddNewVerToPool(group &, vector<loop> &, vector<int> &);
    void _AddNewVer4ToPool(group &, vector<loop> &, vector<int> &);
};

}; // namespace diag

#endif
