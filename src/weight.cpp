#include <iostream>
#include <sstream>
#include <string>
#include <array>
#include "utility/vector.h"
#include "weight.h"

using namespace diag;
using namespace std;

void weight::ReadDiagrams(string FilePrefix)
{
    int count = 0;
    Pool.GPool.clear();
    Pool.VerPool.clear();
    Pool.Ver4Pool.clear();

    while (true)
    {
        count++;
        ifstream DiagFile(FilePrefix + to_string(count) + ".txt");
        if (DiagFile)
        {
            // group Group;
            LOG_INFO("Find " << FilePrefix + to_string(count) + ".txt\n");
            // vector<green> GList;
            istream& DiagFileStream=DiagFile;
            group Group=ReadOneGroup(DiagFileStream, Pool);
            // ReadOneGroup(DiagFile, Group, GList);
        }
        else
            break;
    }
    LOG_INFO("Find " << Pool.GPool.size() << " indepdent green's function.");
    LOG_INFO("Find " << Pool.VerPool.size() << " indepdent interactions.");
    LOG_INFO("Find " << Pool.Ver4Pool.size() << " indepdent 4-vertex.");
}

// vector<green> _ConstructGPool(vector<int> &Permutation, vector<loop> &LoopBasis, vector<int> &GType)
// {
//     vector<green> GPool;
//     int GNum = Permutation.size();
//     for (int i = 0; i < GNum; i++)
//     {
//         cout << "i: " << i << "->" << Permutation[i] << endl;
//         //check the green function from i to Permutation[i]
//         vector<green *> GPool_Type;
//         // select all green in GPooll have the type as GType[i]
//         for (auto &g : GPool)
//             if (g.Type == GType[i])
//                 GPool_Type.push_back(&g);
//         // remove all green's function with different TauIn and TauOut
//         /////   ASSUME Tau(2*i)==Tau(2*i+1) !!!
//         vector<green *> GPool_Tau;
//         for (auto g : GPool_Type)
//         {
//             // cout << g->Type << ":" << g->TauBasis[0] << "->" << g->TauBasis[1] << endl;
//             if (((*g).TauBasis[0] / 2 == i / 2) && ((*g).TauBasis[1] / 2 == Permutation[i] / 2))
//                 GPool_Tau.push_back(g);
//         }

//         // remove all green's function with different LoopBasis
//         vector<green *> GPool_Basis;
//         for (auto g : GPool_Tau)
//         {
//             cout << g->LoopBasis[0] << "," << g->LoopBasis[1] << "," << g->LoopBasis[2] << "," << g->LoopBasis[3] << endl;
//             cout << g->LoopBasis[4] << "," << g->LoopBasis[5] << "," << g->LoopBasis[6] << "," << g->LoopBasis[7] << endl;
//             cout << LoopBasis[i][0] << "," << LoopBasis[i][1] << "," << LoopBasis[i][2] << "," << LoopBasis[i][3] << endl;
//             cout << LoopBasis[i][4] << "," << LoopBasis[i][5] << "," << LoopBasis[i][6] << "," << LoopBasis[i][7] << endl;

//             bool Flag = true;
//             for (int j = 0; j < MaxLoopNum; j++)
//             {
//                 cout << "Compare " << g->LoopBasis[j] << ": " << LoopBasis[i][j] << endl;
//                 if (g->LoopBasis[j] != LoopBasis[i][j])
//                     Flag = false;
//             }
//             if (Flag)
//                 GPool_Basis.push_back(g);
//         }
//         cout << "Basis pool: " << GPool_Basis.size() << "," << MaxLoopNum << endl;
//         //if GPool_Filter is not empty, means the green's function already exists
//         if (GPool_Basis.size() > 0)
//         {

//             cout << "G exits in Pool" << endl;
//             cout << GPool_Basis.size() << endl;
//             ASSERT_ALLWAYS(GPool_Basis.size() == 1, "There are more than two same G in the GPool!");
//             Diagram.GIndex[i] = GPool_Basis[0];
//         }
//         else
//         {
//             cout << "Add G to Pool" << endl;
//             GPool.push_back(MakeG(GType[i], LoopBasis[i], i, Permutation[i]));
//             cout << "Type: " << GPool.back().Type << GPool.back().TauBasis[0] << "->" << GPool.back().TauBasis[1] << endl;
//             green *g = &GPool.back();
//             cout << g->LoopBasis[0] << "," << g->LoopBasis[1] << "," << g->LoopBasis[2] << "," << g->LoopBasis[3] << endl;
//             cout << g->LoopBasis[4] << "," << g->LoopBasis[5] << "," << g->LoopBasis[6] << "," << g->LoopBasis[7] << endl;
//             Diagram.GIndex[i] = &GPool.back();
//         }
//     }
// }

void _AddNewVerToPool(group &Group, vector<vector<int>> &LoopBasisVer, vector<int> &VerType)
{
}
