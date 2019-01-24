#ifndef diagram_H
#define diagram_H

#include <vector>
#include <string>
#include <iostream>
#include "global.h"
#include "utility/utility.h"
#include <array>

namespace diag
{
const size_t MaxBranchNum = 1 << (MaxOrder - 1); //2**(MaxOrder-1)

//column-major two dimensional array
template <class T, size_t ROW, size_t COL>
using matrix = std::array<std::array<T, COL>, ROW>;

typedef std::array<int, MaxLoopNum> loop; //array to store the loop basis for a propagator or interaction line
typedef std::array<double, 2> tau;        //array to store the tau basis (In and Out) for a propagator or interaction line

// G pool to store all basis for Green's functions
struct green
{
    long long int Version; //keep track of the version
    int Type;              //type of each green's function
    loop LoopBasis;        //loop basis for momentum
    tau TauBasis;          //tau basis
    double Weight;         //weight of each green's function
};

struct vertex
{
    // Ver pool to store all basis for interaction lines
    // There two elements, one for direct interaction, another for exchange interaction
    long long int Version;    //keep track of the version
    array<int, 2> Type;       //type of each vertex function
    array<loop, 2> LoopBasis; //loop basis for momentum transfer
    // array<tau, 2> TauBasis;   //tau basis (in and out)
    tau TauBasis;             //tau basis, In and Out
    array<double, 2> Weight;  //weight of each green's function
};

// Ver pool to store all basis for 4-vertex functions
struct vertex4
{
    long long int Version;      //keep track of the version
    int Type;                   //type of each vertex function
    array<green *, 4> Ver4Legs; //the GIndex of four legs of every indepdent 4-vertex
    array<loop, 4> LoopBasis;   //loop basis for momentum transfer
    // tau TauBasis;               //tau basis, Left and Right
    double Weight;              //weight of each green's function
};

struct pool{
    vector<green> GPool;
    vector<vertex> VerPool;
    vector<vertex4> Ver4Pool;
};

struct diagram
{
    double SymFactor;                         //the symmetry factor of a diagram
    array<int, MaxBranchNum> SpinFactor;      //the spin factor of a diagram
    array<green *, 2 * MaxOrder> GIndex;      //the index of all indepdent G
    array<vertex *, 2 * MaxOrder> VerIndex;   //the index of all indepdent interaction
    array<vertex4 *, 2 * MaxOrder> Ver4Index; //the index of all indepdent 4-vertex
};

// A group can be diagrams with different orders,
// or diagrams with same order but have little sign cancellation
// store G and Ver indexes pointing to the corresponding pool
struct group
{
    int ID;
    int HugenNum;             //Number of Hugenholz diagrams in each group
    int Order;                //diagram order of the group
    int LoopNum;              //dimension of loop basis
    int TauNum;               //dimension of tau basis
    int Ver4Num;              //number of 4-vertex
    int GNum;                 //number of G
    vector<diagram> DiagList; //diagrams
};

group ReadOneGroup(istream&, pool&);
}; // namespace diag

#endif