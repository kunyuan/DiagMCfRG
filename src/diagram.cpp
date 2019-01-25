#include "diagram.h"
#include "utility/utility.h"
#include <sstream>

using namespace diag;
using namespace std;

string _GetOneLine(istream &file) {
  // a line only contains space will be skipped
  string buff;
  do {
    getline(file, buff);
  } while (buff.find_first_not_of(' ') == buff.npos);
  return buff;
}

template <typename T> vector<T> _ExtractOneLine(istream &file) {
  // This function extracts the first T type number in one line in the file
  // stream.
  string line;
  getline(file, line);
  stringstream ss;
  vector<T> IntList;
  /* Storing the whole string into string stream */
  ss << line;
  /* Running loop till the end of the stream */
  string temp;
  T found;
  while (!ss.eof()) {
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

vector<loop> _Transpose(vector<vector<int>> &Basis) {
  vector<loop> NewBasis;
  if (Basis.size() == 0)
    return NewBasis;

  for (int i = 0; i < Basis[0].size(); i++) {
    loop TempLoopBasis;
    TempLoopBasis.fill(0);

    for (int j = 0; j < Basis.size(); j++)
      TempLoopBasis[j] = Basis[j][i];

    NewBasis.push_back(TempLoopBasis);
  }
  return NewBasis;
}

green *_AddOneGToPool(vector<green> &GPool, green &Green) {
  vector<green *> GPool_Type;
  // select all green in GPooll have the type as GType[i]
  for (auto &g : GPool)
    if (g.Type == Green.Type)
      GPool_Type.push_back(&g);
  // select all green's function with same TauIn and TauOut
  /////   ASSUME Tau(2*i)==Tau(2*i+1) !!!
  vector<green *> GPool_Tau;
  for (auto g : GPool_Type) {
    // cout << g->Type << ":" << g->TauBasis[0] << "->" << g->TauBasis[1] <<
    // endl;
    // TODO: ugly logic flow
    if (((*g).TauBasis[0] / 2 == Green.TauBasis[0] / 2) &&
        ((*g).TauBasis[1] / 2 == Green.TauBasis[1] / 2) &&
        Green.TauBasis[0] / 2 != 0 && Green.TauBasis[1] / 2 != 0) {
      GPool_Tau.push_back(g);
    }
    if (((*g).TauBasis[0] == Green.TauBasis[0]) &&
        ((*g).TauBasis[1] == Green.TauBasis[1]) && Green.TauBasis[0] / 2 == 0 &&
        Green.TauBasis[1] / 2 == 0) {
      GPool_Tau.push_back(g);
    }
    if (((*g).TauBasis[0] / 2 == Green.TauBasis[0] / 2) &&
        ((*g).TauBasis[1] == Green.TauBasis[1]) && Green.TauBasis[0] / 2 != 0 &&
        Green.TauBasis[1] / 2 == 0) {
      GPool_Tau.push_back(g);
    }
    if (((*g).TauBasis[0] == Green.TauBasis[0]) &&
        ((*g).TauBasis[1] / 2 == Green.TauBasis[1] / 2) &&
        Green.TauBasis[0] / 2 == 0 && Green.TauBasis[1] / 2 != 0) {
      GPool_Tau.push_back(g);
    }
  }
  // select all green's function with same LoopBasis
  vector<green *> GPool_Basis;
  for (auto g : GPool_Tau)
    if (Equal<int>(g->LoopBasis.data(), Green.LoopBasis.data(), MaxLoopNum))
      GPool_Basis.push_back(g);
  // cout << "Basis pool: " << GPool_Basis.size() << "," << MaxLoopNum << endl;
  // if GPool_Filter is not empty, means the green's function already exists
  if (GPool_Basis.size() > 0) {
    ASSERT_ALLWAYS(GPool_Basis.size() == 1,
                   "There are more than two same G in the GPool!");
    return GPool_Basis[0];
  } else {
    GPool.push_back(Green);
    cout << "First: " << ToString(GPool[0]) << endl;
    cout << "Last: " << ToString(GPool.back()) << endl;
    cout << "FirstPointer: " << &GPool[0] << endl;
    cout << "LastPointer: " << &GPool.back() << endl;
    return &(GPool.back());
  }
}

vertex *_AddOneVerToPool(vector<vertex> &VerPool, vertex &Vertex) {
  vector<vertex *> VerPool_Type;
  // select all vertex in VerPooll have the type as Vertex
  for (auto &v : VerPool)
    if (v.Type[0] == Vertex.Type[0])
      if (v.Type[1] == Vertex.Type[1])
        VerPool_Type.push_back(&v);

  vector<vertex *> VerPool_Basis;
  for (auto v : VerPool_Type)
    if (Equal<int>(v->LoopBasis[0].data(), Vertex.LoopBasis[0].data(),
                   MaxLoopNum))
      if (Equal<int>(v->LoopBasis[1].data(), Vertex.LoopBasis[1].data(),
                     MaxLoopNum))
        VerPool_Basis.push_back(v);
  // cout << "Basis pool: " << GPool_Basis.size() << "," << MaxLoopNum << endl;
  // if GPool_Filter is not empty, means the green's function already exists
  if (VerPool_Basis.size() > 0) {
    ASSERT_ALLWAYS(VerPool_Basis.size() == 1,
                   "There are more than two same Vertex in the VerPool!");
    return VerPool_Basis[0];
  } else {
    VerPool.push_back(Vertex);
    return &VerPool.back();
  }
}

vertex4 *_AddOneVer4ToPool(vector<vertex4> &Ver4Pool, vertex4 &Vertex4) {
  vector<vertex4 *> Ver4Pool_Type;
  // select all vertex in VerPooll have the type as Vertex
  for (auto &v : Ver4Pool)
    if (v.Type == Vertex4.Type)
      Ver4Pool_Type.push_back(&v);

  vector<vertex4 *> Ver4Pool_Basis;
  for (auto v : Ver4Pool_Type) {
    bool Flag = true;
    for (int leg = 0; leg < 4; leg++)
      if (!Equal<int>(v->LoopBasis[leg].data(), Vertex4.LoopBasis[leg].data(),
                      MaxLoopNum))
        Flag = false;
    if (Flag)
      Ver4Pool_Basis.push_back(v);
  }
  // cout << "Basis pool: " << GPool_Basis.size() << "," << MaxLoopNum << endl;
  // if GPool_Filter is not empty, means the green's function already exists
  if (Ver4Pool_Basis.size() > 0) {
    ASSERT_ALLWAYS(Ver4Pool_Basis.size() == 1,
                   "There are more than two same vertex4 in the Ver4Pool!");
    return Ver4Pool_Basis[0];
  } else {
    Ver4Pool.push_back(Vertex4);
    return &Ver4Pool.back();
  }
}

vector<green *> _AddAllGToPool(vector<green> &GPool, vector<int> &Permutation,
                               vector<loop> &LoopBasis, vector<int> &GType,
                               int GNum) {
  vector<green *> GIndex;
  for (int i = 0; i < GNum; i++) {
    // construct a new green's function
    green Green;
    Green.Excited = false;
    Green.Version = -1;
    Green.Type = GType[i];
    Green.LoopBasis.fill(0);
    std::copy(LoopBasis[i].begin(), LoopBasis[i].end(),
              Green.LoopBasis.begin());
    Green.TauBasis[IN] = i;
    Green.TauBasis[OUT] = Permutation[i];
    Green.Weight = 1.0e-10;
    // GList.push_back(Green);
    GIndex.push_back(_AddOneGToPool(GPool, Green));
  }
  cout << "All First: " << ToString(GPool[0]) << endl;
  cout << "All Last: " << ToString(GPool.back()) << endl;
  cout << "All First: " << ToString(*(GIndex[0])) << endl;
  cout << "All Last: " << ToString(*(GIndex.back())) << endl;
  cout << GIndex[0] << " vs " << &(GPool[0]) << endl;
  return GIndex;
}

vector<vertex *> _AddAllVerToPool(vector<vertex> &VerPool,
                                  vector<int> &Permutation,
                                  vector<loop> &LoopBasisVer,
                                  vector<int> &VerType, int Ver4Num) {
  vector<vertex *> VerIndex;
  for (int i = 0; i < Ver4Num; i++) {
    int Inidx = 2 * i, Outidx = 2 * i + 1;
    // construct a new green's function
    vertex Vertex;
    Vertex.Excited = {false, false};
    Vertex.Version = -1;
    Vertex.Type = {VerType[Inidx], VerType[Outidx]};
    Vertex.LoopBasis[IN].fill(0);
    Vertex.LoopBasis[OUT].fill(0);
    std::copy(LoopBasisVer[Inidx].begin(), LoopBasisVer[Inidx].end(),
              Vertex.LoopBasis[0].begin());
    std::copy(LoopBasisVer[Outidx].begin(), LoopBasisVer[Outidx].end(),
              Vertex.LoopBasis[0].begin());
    Vertex.Weight = {1.0e-10, 1.0e-10};
    VerIndex.push_back(_AddOneVerToPool(VerPool, Vertex));
  }
  return VerIndex;
}

vector<vertex4 *> _AddAllVer4ToPool(vector<vertex4> &Ver4Pool,
                                    vector<int> &Permutation,
                                    vector<int> &Ver4Legs,
                                    vector<loop> &LoopBasis,
                                    vector<int> &VerType, int Ver4Num) {
  vector<vertex4 *> Ver4Index;
  for (int i = 0; i < Ver4Num; i++) {
    int Inidx = 2 * i, Outidx = 2 * i + 1;
    // construct a new green's function
    vertex4 Vertex4;
    Vertex4.Excited = false;
    Vertex4.Version = 0;
    Vertex4.Type = VerType[Inidx];
    for (int leg = 0; leg < 4; leg++) {
      Vertex4.LoopBasis[leg].fill(0);
      int gidx = Ver4Legs[leg];
      std::copy(LoopBasis[gidx].begin(), LoopBasis[gidx].end(),
                Vertex4.LoopBasis[leg].begin());
    }
    Vertex4.Weight = 1.0e-10;
    Ver4Index.push_back(_AddOneVer4ToPool(Ver4Pool, Vertex4));
  }
  return Ver4Index;
}

diagram ReadOneDiagram(istream &DiagFile, pool &Pool, int Order, int LoopNum,
                       int GNum, int Ver4Num) {
  string buff;
  diagram Diagram;

  Diagram.GIndex.fill(nullptr);
  Diagram.VerIndex.fill(nullptr);
  Diagram.Ver4Index.fill(nullptr);

  //////// Diagram Topology  ////////////////////////
  buff = _GetOneLine(DiagFile); // title
  auto Permutation = _ExtractOneLine<int>(DiagFile);

  //////// Propagator type //////////////////
  buff = _GetOneLine(DiagFile); // title
  auto GType = _ExtractOneLine<int>(DiagFile);

  //////// symmetry factor //////////////////
  buff = _GetOneLine(DiagFile); // title
  Diagram.SymFactor = _ExtractOneLine<double>(DiagFile)[0];

  /////// Loop Basis  /////////////////////////
  buff = _GetOneLine(DiagFile); // title
  vector<vector<int>> TransposedLoopBasis;
  for (int j = 0; j < LoopNum; j++)
    TransposedLoopBasis.push_back(_ExtractOneLine<int>(DiagFile));
  vector<loop> LoopBasis = _Transpose(TransposedLoopBasis);

  /////// 4 legs of 4-ver  /////////////////////////
  buff = _GetOneLine(DiagFile); // title
  vector<int> Ver4Legs;
  if (Ver4Num > 0)
    Ver4Legs = _ExtractOneLine<int>(DiagFile);

  /////// Interaction Type  /////////////////////////
  buff = _GetOneLine(DiagFile); // title
  vector<int> VerType;
  if (Ver4Num > 0)
    VerType = _ExtractOneLine<int>(DiagFile);

  /////// Interaction Loop Basis  ////////////////////
  buff = _GetOneLine(DiagFile); // title
  vector<vector<int>> TransposedLoopBasisVer;
  if (Ver4Num > 0)
    for (int j = 0; j < LoopNum; j++)
      TransposedLoopBasisVer.push_back(_ExtractOneLine<int>(DiagFile));
  vector<loop> LoopBasisVer = _Transpose(TransposedLoopBasisVer);

  /////// Spin Factor  ////////////////////
  buff = _GetOneLine(DiagFile); // title
  auto SpinFactor = _ExtractOneLine<double>(DiagFile);
  copy(SpinFactor.begin(), SpinFactor.end(),
       Diagram.SpinFactor.begin()); // copy spin factor into the group member

  ////////   Add G to GPool ////////////////////////
  vector<green *> GIndex =
      _AddAllGToPool(Pool.GPool, Permutation, LoopBasis, GType, GNum);
  ASSERT_ALLWAYS(GIndex.size() == GNum,
                 "Number of Green's function does not match!");
  copy(GIndex.begin(), GIndex.end(), Diagram.GIndex.begin());

  cout << "in diagram" << endl;
  cout << ToString(*(Diagram.GIndex[0])) << endl;
  cout << Diagram.GIndex[0] << endl;

  if (Ver4Num > 0) {
    //////  Add 4-Vertex to Ver4Pool /////////
    vector<vertex4 *> Ver4Index = _AddAllVer4ToPool(
        Pool.Ver4Pool, Permutation, Ver4Legs, LoopBasis, VerType, Ver4Num);
    copy(Ver4Index.begin(), Ver4Index.end(), Diagram.Ver4Index.begin());
    //////  Add Vertex to VerPool /////////
    vector<vertex *> VerIndex = _AddAllVerToPool(
        Pool.VerPool, Permutation, LoopBasisVer, VerType, Ver4Num);
    copy(VerIndex.begin(), VerIndex.end(), Diagram.VerIndex.begin());
  }

  return Diagram;
}

group diag::ReadOneGroup(istream &DiagFile, pool &Pool) {
  group Group;
  Group.HugenNum = _ExtractOneLine<int>(DiagFile)[0];
  Group.Order = _ExtractOneLine<int>(DiagFile)[0];
  string buff = _GetOneLine(DiagFile); // skip the line for the group type

  Group.LoopNum = Group.Order + 1;
  Group.InternalLoopNum = Group.Order;
  Group.TauNum = Group.Order * 2;
  Group.InternalTauNum = Group.Order * 2 - 2;
  Group.GNum = Group.Order * 2;
  Group.Ver4Num = Group.Order - 1;

  for (int i = 0; i < Group.HugenNum; i++) {
    diagram Diagram = ReadOneDiagram(DiagFile, Pool, Group.Order, Group.LoopNum,
                                     Group.GNum, Group.Ver4Num);
    Diagram.ID = i;
    Group.DiagList.push_back(Diagram);
    // for (int i = 0; i < Group.GNum; i++)
    //   cout << ToString(*(Group.DiagList.back().GIndex[i])) << endl;
  }
  diag::Test(Group);
  //   cout << "Group" << endl;
  //   for (int i = 0; i < Group.GNum; i++)
  cout << "in Group" << endl;
  cout << ToString(Group) << endl;
  cout << ToString(*(Group.DiagList[0].GIndex[0])) << endl;

  return Group;
}

std::string ToString(const diagram &Diag) {
  std::ostringstream oss;
  oss << "DiagID: " << Diag.ID << "\n Weight=" << Diag.Weight
      << "\n SymFactor=" << Diag.SymFactor << "\n SpinFactor="
      << ToString<double>(Diag.SpinFactor.data(),
                          (size_t)Diag.SpinFactor.size())
      << endl;
  return oss.str();
}
std::string ToString(const group &Group) {
  std::ostringstream oss;
  oss << "GroupID: " << Group.ID << "\n Weight=" << Group.Weight
      << "\n HugenNum=" << Group.HugenNum << endl;
  return oss.str();
};
std::string ToString(const green &G) {
  std::ostringstream oss;
  oss << "Version: " << G.Version << "\n Type=" << G.Type
      << "\n Weight=" << G.Weight
      << "\n TauBasis=" << ToString<int>(G.TauBasis.data(), 2)
      << "\n LoopBasis="
      << ToString<int>(G.LoopBasis.data(), (size_t)G.LoopBasis.size()) << endl;
  return oss.str();
};
std::string ToString(const vertex &);
std::string ToString(const vertex4 &);

void diag::Test(group &Group) {
  for (auto &d : Group.DiagList) {
    for (auto i = 0; i < Group.GNum; i++) {
      CHECKNULL(d.GIndex[i]);
    }
    for (auto i = 0; i < Group.Ver4Num; i++) {
      CHECKNULL(d.VerIndex[i]);
      CHECKNULL(d.Ver4Index[i]);
    }
  }
}
