#include <iostream>
#include <sstream>
#include <string>
#include <array>
#include "utility/vector.h"
#include "weight.h"

using namespace diag;
using namespace std;

template <typename T>
vector<T> ExtractOneLine(ifstream &file)
{
    //This function extracts the first T type number in one line in the file stream.
    string line;
    getline(file, line);
    stringstream ss;
    vector<T> IntList;
    /* Storing the whole string into string stream */
    ss << line;
    /* Running loop till the end of the stream */
    string temp;
    T found;
    while (!ss.eof())
    {
        /* extracting word by word from stream */
        ss >> temp;
        /* Checking the given word is integer or not */
        if (stringstream(temp) >> found)
            IntList.push_back(found);
        /* To save from space at the end of string */
        temp = "";
    }
    return IntList;
}

void weight::ReadDiagrams(string FilePrefix)
{
    int count = 0;
    while (true)
    {
        count++;
        ifstream DiagFile(FilePrefix + to_string(count) + ".txt");
        if (DiagFile)
        {
            group Group;
            LOG_INFO("Find " << FilePrefix + to_string(count) + ".txt\n");
            _ReadOneGroup(DiagFile, Group);
        }
        else
            break;
    }
}

void weight::_ReadOneGroup(ifstream &DiagFile, group &Group)
{
    string buff;
    Group.HugenNum = ExtractOneLine<int>(DiagFile)[0];
    Group.Order = ExtractOneLine<int>(DiagFile)[0];
    Group.LoopNum = Group.Order + 1;
    Group.TauNum = Group.Order * 2;
    Group.GNum = Group.Order * 2;
    Group.Ver4Num = Group.Order - 1;

    getline(DiagFile, buff); //skip the line for the diagram type
    getline(DiagFile, buff); //skip the blank line

    for (int i = 0; i < Group.HugenNum; i++)
    {
        //////// Diagram Topology  ////////////////////////
        getline(DiagFile, buff); //title
        auto Permutation = ExtractOneLine<int>(DiagFile);

        //////// Propagator type //////////////////
        getline(DiagFile, buff); //title
        auto GType = ExtractOneLine<int>(DiagFile);

        //////// symmetry factor //////////////////
        getline(DiagFile, buff); //title
        Group.SymFactor[i] = ExtractOneLine<double>(DiagFile)[0];

        /////// Loop Basis  /////////////////////////
        getline(DiagFile, buff); //title
        vector<vector<int>> LoopBasis;
        for (int j = 0; j < Group.LoopNum; j++)
            LoopBasis.push_back(ExtractOneLine<int>(DiagFile));

        /////// 4 legs of 4-ver  /////////////////////////
        getline(DiagFile, buff); //title
        vector<int> Ver4Legs;
        if (Group.Order > 1)
            Ver4Legs = ExtractOneLine<int>(DiagFile);

        /////// Interaction Type  /////////////////////////
        getline(DiagFile, buff); //title
        vector<int> VerType;
        if (Group.Order > 1)
            VerType = ExtractOneLine<int>(DiagFile);

        /////// Interaction Loop Basis  ////////////////////
        getline(DiagFile, buff); //title
        vector<vector<int>> LoopBasisVer;
        if (Group.Order > 1)
            for (int j = 0; j < Group.LoopNum; j++)
                LoopBasisVer.push_back(ExtractOneLine<int>(DiagFile));

        /////// Spin Factor  ////////////////////
        getline(DiagFile, buff); //title
        auto SpinFactor = ExtractOneLine<int>(DiagFile);
        // AssignFromTo(SpinFactor.data(), Group.SpinFactor[i].data(), SpinFactor.size());
        copy(SpinFactor.begin(), SpinFactor.end(), Group.SpinFactor[i].begin()); //copy spin factor into the group member

        getline(DiagFile, buff); //blank

        // _AddNewGToPool(Group, Permutation, LoopBasis, GType);
        // _AddNewVerToPool(Group, LoopBasisVer, VerType);
        // _AddNewVer4ToPool(Group, LoopBasis, VerType);
    }
    return;
}

void _AddNewGToPool(group &Group, vector<int> &Permutation, vector<vector<int>> &LoopBasis, vector<int> &GType)
{
    // for (int i = 0; i < Group.TauNum; i++)
    // {
    // }
}

void _AddNewVerToPool(group &Group, vector<vector<int>> &LoopBasisVer, vector<int> &VerType)
{
}
