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
    // cout << "size: " << IntList.size() << endl;
    // cout << IntList[0] << endl;
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
    LOG_INFO("Find " << GPool.size() << " indepdent green's function.");
    LOG_INFO("Find " << VerPool.size() << " indepdent interactions.");
    LOG_INFO("Find " << Ver4Pool.size() << " indepdent 4-vertex.");
}

vector<vector<int>> Transpose(vector<vector<int>> &Basis)
{
    vector<vector<int>> NewBasis;
    if (Basis.size() == 0)
        return Basis;
    for (int i = 0; i < Basis[0].size(); i++)
    {
        NewBasis.push_back(vector<int>());
        for (int j = 0; j < Basis.size(); j++)
            NewBasis[i].push_back(Basis[j][i]);
        for (int j = Basis.size(); j < MaxLoopNum; j++)
            NewBasis[i].push_back(0); //pending zero
    }
    return NewBasis;
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
        Group.DiagList.push_back(diagram()); //add a new diagram structure

        //////// Diagram Topology  ////////////////////////
        getline(DiagFile, buff); //title
        auto Permutation = ExtractOneLine<int>(DiagFile);

        //////// Propagator type //////////////////
        getline(DiagFile, buff); //title
        auto GType = ExtractOneLine<int>(DiagFile);
        cout << "GType" << GType[0] << endl;

        //////// symmetry factor //////////////////
        getline(DiagFile, buff); //title
        Group.DiagList[i].SymFactor = ExtractOneLine<double>(DiagFile)[0];

        /////// Loop Basis  /////////////////////////
        getline(DiagFile, buff); //title
        vector<vector<int>> LoopBasis;
        for (int j = 0; j < Group.LoopNum; j++)
            LoopBasis.push_back(ExtractOneLine<int>(DiagFile));
        LoopBasis = Transpose(LoopBasis);

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
        LoopBasisVer = Transpose(LoopBasisVer);

        /////// Spin Factor  ////////////////////
        getline(DiagFile, buff); //title
        auto SpinFactor = ExtractOneLine<int>(DiagFile);
        // AssignFromTo(SpinFactor.data(), Group.SpinFactor[i].data(), SpinFactor.size());
        copy(SpinFactor.begin(), SpinFactor.end(), Group.DiagList[i].SpinFactor.begin()); //copy spin factor into the group member

        getline(DiagFile, buff); //blank

        _AddNewGToPool(Group.DiagList[i], Permutation, LoopBasis, GType);
        // _AddNewVerToPool(Group, LoopBasisVer, VerType);
        // _AddNewVer4ToPool(Group, LoopBasis, VerType);
    }
    return;
}

void MakeG(green &Diagram, int GType, vector<int> &LoopBasis, int TauIn, int TauOut)
{
    // cout << "Make new G" << endl;
    Diagram.Version = 0;
    Diagram.Type = GType;
    Diagram.LoopBasis.fill(0); //all basis will initilize with 0!
    // cout << "copy" << Diagram.LoopBasis.size() << "," << LoopBasis.size() << endl;
    copy(LoopBasis.begin(), LoopBasis.end(), Diagram.LoopBasis.begin());
    // cout << "copy finished" << endl;
    Diagram.TauBasis[IN] = TauIn;
    Diagram.TauBasis[OUT] = TauOut;
    Diagram.Weight = 1.0e-10;

    // cout << "Type: " << GType << ", Tau:" << TauIn << "," << TauOut << endl;
}

void weight::_AddNewGToPool(diagram &Diagram, vector<int> &Permutation, vector<vector<int>> &LoopBasis, vector<int> &GType)
{
    int GNum = Permutation.size();
    for (int i = 0; i < GNum; i++)
    {
        cout << "i: " << i << "->" << Permutation[i] << endl;
        //check the green function from i to Permutation[i]
        vector<green *> GPool_Type;
        // select all green in GPooll have the type as GType[i]
        for (auto &g : GPool)
            if (g.Type == GType[i])
                GPool_Type.push_back(&g);
        // remove all green's function with different TauIn and TauOut
        /////   ASSUME Tau(2*i)==Tau(2*i+1) !!!
        vector<green *> GPool_Tau;
        for (auto g : GPool_Type)
        {
            // cout << g->Type << ":" << g->TauBasis[0] << "->" << g->TauBasis[1] << endl;
            if (((*g).TauBasis[0] / 2 == i / 2) && ((*g).TauBasis[1] / 2 == Permutation[i] / 2))
                GPool_Tau.push_back(g);
        }

        // remove all green's function with different LoopBasis
        vector<green *> GPool_Basis;
        for (auto g : GPool_Tau)
        {
            cout << g->LoopBasis[0] << "," << g->LoopBasis[1] << "," << g->LoopBasis[2] << "," << g->LoopBasis[3] << endl;
            cout << g->LoopBasis[4] << "," << g->LoopBasis[5] << "," << g->LoopBasis[6] << "," << g->LoopBasis[7] << endl;
            cout << LoopBasis[i][0] << "," << LoopBasis[i][1] << "," << LoopBasis[i][2] << "," << LoopBasis[i][3] << endl;
            cout << LoopBasis[i][4] << "," << LoopBasis[i][5] << "," << LoopBasis[i][6] << "," << LoopBasis[i][7] << endl;

            bool Flag = true;
            for (int j = 0; j < MaxLoopNum; j++)
            {
                cout << "Compare " << g->LoopBasis[j] << ": " << LoopBasis[i][j] << endl;
                if (g->LoopBasis[j] != LoopBasis[i][j])
                    Flag = false;
            }
            if (Flag)
                GPool_Basis.push_back(g);
        }
        cout << "Basis pool: " << GPool_Basis.size() << "," << MaxLoopNum << endl;
        //if GPool_Filter is not empty, means the green's function already exists
        if (GPool_Basis.size() > 0)
        {

            cout << "G exits in Pool" << endl;
            cout << GPool_Basis.size() << endl;
            ASSERT_ALLWAYS(GPool_Basis.size() == 1, "There are more than two same G in the GPool!");
            Diagram.GIndex[i] = GPool_Basis[0];
        }
        else
        {
            cout << "Add G to Pool" << endl;
            GPool.push_back(green());
            MakeG(GPool.back(), GType[i], LoopBasis[i], i, Permutation[i]);
            cout << "Type: " << GPool.back().Type << GPool.back().TauBasis[0] << "->" << GPool.back().TauBasis[1] << endl;
            green *g = &GPool.back();
            cout << g->LoopBasis[0] << "," << g->LoopBasis[1] << "," << g->LoopBasis[2] << "," << g->LoopBasis[3] << endl;
            cout << g->LoopBasis[4] << "," << g->LoopBasis[5] << "," << g->LoopBasis[6] << "," << g->LoopBasis[7] << endl;
            Diagram.GIndex[i] = &GPool.back();
        }
    }
}

void _AddNewVerToPool(group &Group, vector<vector<int>> &LoopBasisVer, vector<int> &VerType)
{
}
