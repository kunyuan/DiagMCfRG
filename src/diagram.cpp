#include "diagram.h"

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

vector<loop> Transpose(vector<vector<int>> &Basis)
{
    vector<loop> NewBasis;
    if (Basis.size() == 0)
        return NewBasis;
    for (int i = 0; i < Basis[0].size(); i++)
    {
        loop TempLoopBasis;
        TempLoopBasis.fill(0);
        for (int j = 0; j < Basis.size(); j++)
            TempLoopBasis[j] = Basis[j][i];
        NewBasis.push_back(TempLoopBasis);
    }
    return NewBasis;
}

vector<green> _ConstructGPool(vector<int> &Permutation, vector<loop> &LoopBasis, vector<int> &GType)
{
    vector<green> GPool;
    int GNum = Permutation.size();
    for (int i = 0; i < GNum; i++)
    {
        //construct a new green's function
        green Green;
        Green.Version = 0;
        Green.Type = GType[i];
        Green.LoopBasis.fill(0);
        copy(LoopBasis[i].begin(), LoopBasis[i].end(), Green.LoopBasis.begin());
        Green.TauBasis[IN] = i;
        Green.TauBasis[OUT] = Permutation[i];
        Green.Weight = 1.0e-10;
        GPool.push_back(Green);
    }
    return GPool;
}

diagram ReadOneDiagram(ifstream &DiagFile, int Order, int LoopNum)
{
    string buff;
    //////// Diagram Topology  ////////////////////////
    getline(DiagFile, buff); //title
    auto Permutation = ExtractOneLine<int>(DiagFile);

    //////// Propagator type //////////////////
    getline(DiagFile, buff); //title
    auto GType = ExtractOneLine<int>(DiagFile);
    cout << "GType" << GType[0] << endl;

    //////// symmetry factor //////////////////
    getline(DiagFile, buff); //title
    auto SymFactor = ExtractOneLine<double>(DiagFile)[0];

    /////// Loop Basis  /////////////////////////
    getline(DiagFile, buff); //title
    vector<vector<int>> TransposedLoopBasis;
    for (int j = 0; j < LoopNum; j++)
        TransposedLoopBasis.push_back(ExtractOneLine<int>(DiagFile));
    vector<loop> LoopBasis = Transpose(TransposedLoopBasis);

    /////// 4 legs of 4-ver  /////////////////////////
    getline(DiagFile, buff); //title
    vector<int> Ver4Legs;
    if (Order > 1)
        Ver4Legs = ExtractOneLine<int>(DiagFile);

    /////// Interaction Type  /////////////////////////
    getline(DiagFile, buff); //title
    vector<int> VerType;
    if (Order > 1)
        VerType = ExtractOneLine<int>(DiagFile);

    /////// Interaction Loop Basis  ////////////////////
    getline(DiagFile, buff); //title
    vector<vector<int>> TransposedLoopBasisVer;
    if (Order > 1)
        for (int j = 0; j < LoopNum; j++)
            TransposedLoopBasisVer.push_back(ExtractOneLine<int>(DiagFile));
    vector<loop> LoopBasisVer = Transpose(TransposedLoopBasisVer);

    /////// Spin Factor  ////////////////////
    getline(DiagFile, buff); //title
    auto SpinFactor = ExtractOneLine<int>(DiagFile);

    //construct a diagram
    diagram Diagram;
    copy(SpinFactor.begin(), SpinFactor.end(), Diagram.SpinFactor.begin()); //copy spin factor into the group member
    Diagram.SymFactor = SymFactor;
    //TODO: add other members of diagram
    return Diagram;
}

void diag::_ReadOneGroup(ifstream &DiagFile, group &Group, vector<green> &GPool)
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
        Group.DiagList.push_back(ReadOneDiagram(DiagFile, Group.Order, Group.LoopNum)); //add a new diagram structure
        getline(DiagFile, buff);                                                        //blank

        // _AddNewGToPool(Group.DiagList[i], Permutation, LoopBasis, GType);
        // _AddNewVerToPool(Group, LoopBasisVer, VerType);
        // _AddNewVer4ToPool(Group, LoopBasis, VerType);
    }
    return;
}
