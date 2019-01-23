#ifndef weight_H
#define weight_H

#include <vector>
#include <string>
#include <iostream>
#include "global.h"
#include "utility/utility.h"
#include <Eigen/Dense>
extern parameter Para;

namespace diag
{
using namespace Eigen;
const int MaxBranchNum = 1 << (MaxOrder - 1); //2**(MaxOrder-1)
// A group can be diagrams with different orders,
// or diagrams with same order but have little sign cancellation
// store G and Ver indexes pointing to the corresponding pool
struct group
{
    int HugenNum;                             //Number of Hugenholz diagrams in each group
    int Order;                                //diagram order of the group
    int LoopNum;                              //dimension of loop basis
    int TauNum;                               //dimension of tau basis
    int Ver4Num;                              //number of 4-vertex
    int GNum;                                 //number of G
    double SymFactor[MaxDiagNum];             //the symmetry factor of a diagram
    int SpinFactor[MaxDiagNum][MaxBranchNum]; //the spin factor of a diagram
    int GIndex[MaxDiagNum][2 * MaxOrder];     //the index of all indepdent G
    int VerIndex[MaxDiagNum][2 * MaxOrder];   //the index of all indepdent interaction
    int Ver4Index[MaxDiagNum][MaxGroupNum];   //the index of all indepdent 4-vertex
};

// G pool to store all basis for Green's functions
struct gPool
{
    int Type[MaxGPoolSize]; //type of each green's function
    // int LoopBasis[MaxGPoolSize][MaxOrder + 1]; //loop basis for momentum
    Matrix<int, MaxGPoolSize, MaxOrder + 1> LoopBasis; //loop basis for momentum
    double TauBasis[MaxGPoolSize][2];                  //tau basis
    double Weight[MaxGPoolSize];                       //weight of each green's function
};

// Ver pool to store all basis for interaction lines
struct verPool
{
    int Type[2 * MaxVerPoolSize];                    //type of each vertex function
    int LoopBasis[2 * MaxVerPoolSize][MaxOrder + 1]; //loop basis for momentum transfer
    double TauBasis[2 * MaxVerPoolSize][2];          //tau basis, In and Out
    double Weight[2 * MaxVerPoolSize];               //weight of each green's function
};

// Ver pool to store all basis for 4-vertex functions
struct ver4Pool
{
    int Type[MaxVerPoolSize];                          //type of each vertex function
    int Ver4Legs[4][MaxVerPoolSize];                   //the GIndex of four legs of every indepdent 4-vertex
    double Weight[MaxVerPoolSize];                     //weight of each green's function
    int LoopBasis[4 * MaxVerPoolSize][MaxOrder + 1];   //loop basis for momentum transfer
    double TauBasis[4 * MaxVerPoolSize][2 * MaxOrder]; //tau basis, In and Out
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
    std::vector<group> GroupList;
    gPool GPool;
    verPool VerPool;
    ver4Pool Ver4Pool;
    void _ReadOneGroup(ifstream &, group &);
    void _AddNewGToPool(group &, vector<int> &, vector<vector<int>> &, vector<int> &);
    void _AddNewVerToPool(group &, vector<vector<int>> &, vector<int> &);
    void _AddNewVer4ToPool(group &, vector<vector<int>> &, vector<int> &);
};

}; // namespace diag

#endif
