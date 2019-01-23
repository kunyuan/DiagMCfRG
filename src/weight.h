#ifndef weight_H
#define weight_H

#include <vector>
#include <string>
#include <iostream>
#include "global.h"
#include "utility/utility.h"
#include <array>
extern parameter Para;

namespace diag
{
using namespace std;
const size_t MaxBranchNum = 1 << (MaxOrder - 1); //2**(MaxOrder-1)

//column-major two dimensional array
template <class T, size_t ROW, size_t COL>
using matrix = std::array<std::array<T, COL>, ROW>;

// G pool to store all basis for Green's functions
struct green
{
    int Type;                           //type of each green's function
    array<int, MaxOrder + 1> LoopBasis; //loop basis for momentum
    array<double, 2> TauBasis;          //tau basis
    double Weight;                      //weight of each green's function
};

struct vertex
{
    // Ver pool to store all basis for interaction lines
    // There two elements, one for direct interaction, another for exchange interaction
    array<int, 2> Type;                     //type of each vertex function
    matrix<int, 2, MaxOrder + 1> LoopBasis; //loop basis for momentum transfer
    array<double, 2> TauBasis;              //tau basis, In and Out
    array<double, 2> Weight;                //weight of each green's function
};

// Ver pool to store all basis for 4-vertex functions
struct vertex4
{
    int Type;                               //type of each vertex function
    array<green *, 4> Ver4Legs;             //the GIndex of four legs of every indepdent 4-vertex
    matrix<int, 4, MaxOrder + 1> LoopBasis; //loop basis for momentum transfer
    array<double, 2> TauBasis;              //tau basis, Left and Right
    double Weight;                          //weight of each green's function
};

// A group can be diagrams with different orders,
// or diagrams with same order but have little sign cancellation
// store G and Ver indexes pointing to the corresponding pool
struct group
{
    int HugenNum;                                         //Number of Hugenholz diagrams in each group
    int Order;                                            //diagram order of the group
    int LoopNum;                                          //dimension of loop basis
    int TauNum;                                           //dimension of tau basis
    int Ver4Num;                                          //number of 4-vertex
    int GNum;                                             //number of G
    array<double, MaxDiagNum> SymFactor;                  //the symmetry factor of a diagram
    matrix<int, MaxDiagNum, MaxBranchNum> SpinFactor;     //the spin factor of a diagram
    matrix<green *, MaxDiagNum, 2 * MaxOrder> GIndex;     //the index of all indepdent G
    matrix<vertex *, MaxDiagNum, 2 * MaxOrder> VerIndex;  //the index of all indepdent interaction
    matrix<vertex4 *, MaxDiagNum, MaxGroupNum> Ver4Index; //the index of all indepdent 4-vertex
};


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
    void _ReadOneGroup(ifstream &, group &);
    void _AddNewGToPool(group &, vector<int> &, vector<vector<int>> &, vector<int> &);
    void _AddNewVerToPool(group &, vector<vector<int>> &, vector<int> &);
    void _AddNewVer4ToPool(group &, vector<vector<int>> &, vector<int> &);
};

}; // namespace diag

#endif
