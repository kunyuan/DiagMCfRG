#include "diagram.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/utility.h"
#include <iostream>
#include <sstream>

using namespace diag;
using namespace std;

string _GetOneLine(istream &file) {
  // a line only contains space will be skipped
  string buff;
  do {
    getline(file, buff);
    if (file.bad())
      ABORT("Fail to read the file!");
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

vector<loop> _Transpose(vector<vector<double>> &Basis) {
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

green *_AddOneGToPool(pool &Pool, green &Green) {
  vector<green *> GPool_Type;
  // select all green in GPooll have the type as GType[i]
  for (int i = 0; i < Pool.GPoolSize; i++) {
    green *g = &(Pool.GPool[i]);
    if (g->Type == Green.Type)
      GPool_Type.push_back(g);
  }
  // select all green's function with same TauIn and TauOut
  vector<green *> GPool_Tau;
  for (auto g : GPool_Type) {
    if (g->TauBasis[IN] == Green.TauBasis[IN] &&
        g->TauBasis[OUT] == Green.TauBasis[OUT]) {
      GPool_Tau.push_back(g);
    }
  }
  // select all green's function with same LoopBasis
  vector<green *> GPool_Basis;
  for (auto g : GPool_Tau)
    if (Equal<double>(g->LoopBasis.data(), Green.LoopBasis.data(), MaxLoopNum))
      GPool_Basis.push_back(g);

  //   cout << "Basis pool: " << GPool_Basis.size() << "," << MaxLoopNum <<
  //   endl;

  // if GPool_Filter is not empty, means the green's function already exists
  if (GPool_Basis.size() > 0) {
    for (auto &g : GPool_Basis) {
      //   cout << ToString(Green) << endl;
      //   cout << ToString(*g) << endl;
    }
    ASSERT_ALLWAYS(GPool_Basis.size() == 1,
                   "There are more than two same G in the GPool!");
    return GPool_Basis[0];
  } else {
    ASSERT_ALLWAYS(Pool.GPoolSize < MaxGPoolSize, "MaxGPoolSize is too small!");
    Pool.GPoolSize++;
    Pool.GPool[Pool.GPoolSize - 1] = Green;
    return &Pool.GPool[Pool.GPoolSize - 1];
  }
}

vertex *_AddOneVerToPool(pool &Pool, vertex &Vertex) {
  vector<vertex *> VerPool_Type;
  // select all vertex in VerPooll have the type as Vertex

  for (int i = 0; i < Pool.VerPoolSize; i++) {
    vertex *v = &(Pool.VerPool[i]);
    if (v->Type[0] == Vertex.Type[0])
      if (v->Type[1] == Vertex.Type[1])
        VerPool_Type.push_back(v);
  }

  vector<vertex *> VerPool_Basis;
  for (auto v : VerPool_Type)
    if (Equal<double>(v->LoopBasis[0].data(), Vertex.LoopBasis[0].data(),
                      MaxLoopNum))
      if (Equal<double>(v->LoopBasis[1].data(), Vertex.LoopBasis[1].data(),
                        MaxLoopNum))
        VerPool_Basis.push_back(v);
  // cout << "Basis pool: " << GPool_Basis.size() << "," << MaxLoopNum << endl;
  // if GPool_Filter is not empty, means the green's function already exists
  if (VerPool_Basis.size() > 0) {
    ASSERT_ALLWAYS(VerPool_Basis.size() == 1,
                   "There are more than two same Vertex in the VerPool!");
    return VerPool_Basis[0];
  } else {
    // add new vertex4 to pool
    ASSERT_ALLWAYS(Pool.VerPoolSize < MaxVerPoolSize,
                   "MaxVerPoolSize is too small!");
    Pool.VerPoolSize++;
    Pool.VerPool[Pool.VerPoolSize - 1] = Vertex;
    return &Pool.VerPool[Pool.VerPoolSize - 1];
  }
}

vertex4 *_AddOneInteractionToPool(pool &Pool, vertex4 &Vertex4) {
  vector<vertex4 *> Ver4Pool_Type;
  // select all vertex in VerPooll have the type as Vertex

  for (int i = 0; i < Pool.Ver4PoolSize; i++) {
    vertex4 *v = &(Pool.Ver4Pool[i]);
    if (v->Type[0] == Vertex4.Type[0])
      if (v->Type[1] == Vertex4.Type[1])
        Ver4Pool_Type.push_back(v);
  }

  vector<vertex4 *> Ver4Pool_Basis;
  for (auto v : Ver4Pool_Type)
    if (Equal<double>(v->IntLoopBasis[0].data(), Vertex4.IntLoopBasis[0].data(),
                      MaxLoopNum) &&
        Equal<double>(v->IntLoopBasis[1].data(), Vertex4.IntLoopBasis[1].data(),
                      MaxLoopNum))
      Ver4Pool_Basis.push_back(v);
  // cout << "Basis pool: " << GPool_Basis.size() << "," << MaxLoopNum << endl;
  // if GPool_Filter is not empty, means the green's function already exists
  if (Ver4Pool_Basis.size() > 0) {
    ASSERT_ALLWAYS(Ver4Pool_Basis.size() == 1,
                   "There are more than two same Vertex in the VerPool!");
    return Ver4Pool_Basis[0];
  } else {
    // add new vertex4 to pool
    ASSERT_ALLWAYS(Pool.Ver4PoolSize < MaxVerPoolSize,
                   "MaxVerPoolSize is too small!");
    Pool.Ver4PoolSize++;
    Pool.Ver4Pool[Pool.Ver4PoolSize - 1] = Vertex4;
    return &Pool.Ver4Pool[Pool.Ver4PoolSize - 1];
  }
}

vertex4 *_AddOneVer4ToPool(pool &Pool, vertex4 &Vertex4) {
  vector<vertex4 *> Ver4Pool_Type;
  // select all vertex in VerPooll have the type as Vertex
  for (int i = 0; i < Pool.Ver4PoolSize; i++) {
    vertex4 *v = &(Pool.Ver4Pool[i]);
    if (v->Type[0] == Vertex4.Type[0])
      if (v->Type[1] == Vertex4.Type[1])
        Ver4Pool_Type.push_back(v);
  }

  vector<vertex4 *> Ver4Pool_Basis;
  for (auto v : Ver4Pool_Type) {
    bool Flag = true;
    for (int leg = 0; leg < 4; leg++)
      if (!Equal<double>(v->LoopBasis[leg].data(),
                         Vertex4.LoopBasis[leg].data(), MaxLoopNum))
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
    // add new vertex4 to pool
    ASSERT_ALLWAYS(Pool.Ver4PoolSize < MaxVerPoolSize,
                   "MaxVerPoolSize is too small!");
    Pool.Ver4PoolSize++;
    Pool.Ver4Pool[Pool.Ver4PoolSize - 1] = Vertex4;
    return &Pool.Ver4Pool[Pool.Ver4PoolSize - 1];
  }
}

vector<green *> _AddAllGToPool(pool &Pool, vector<tau> &VerBasis,
                               vector<loop> &LoopBasis, vector<int> &GType,
                               int GNum) {
  vector<green *> GIndex;
  for (int i = 0; i < GNum; i++) {
    // construct a new green's function
    green Green;
    Green.Type = GType[i];
    Green.LoopBasis.fill(0);
    std::copy(LoopBasis[i].begin(), LoopBasis[i].end(),
              Green.LoopBasis.begin());
    Green.TauBasis = VerBasis[i];
    // GList.push_back(Green);
    GIndex.push_back(_AddOneGToPool(Pool, Green));
  }
  //   cout << "All First: " << ToString(GPool[0]) << endl;
  //   cout << "All Last: " << ToString(GPool.back()) << endl;
  //   cout << "All First: " << ToString(*(GIndex[0])) << endl;
  //   cout << "All Last: " << ToString(*(GIndex.back())) << endl;
  //   cout << GIndex[0] << " vs " << &(GPool[0]) << endl;
  return GIndex;
}

vector<vertex *> _AddAllVerToPool(pool &Pool, vector<tau> &VerBasis,
                                  vector<loop> &LoopBasisVer,
                                  vector<int> &VerType, int Ver4Num) {
  vector<vertex *> VerIndex;
  for (int i = 0; i < Ver4Num; i++) {
    int Inidx = 2 * i, Outidx = 2 * i + 1;
    // construct a new vertex function
    vertex Vertex;
    Vertex.Type = {VerType[Inidx], VerType[Outidx]};
    Vertex.LoopBasis[IN].fill(0);
    Vertex.LoopBasis[OUT].fill(0);
    std::copy(LoopBasisVer[Inidx].begin(), LoopBasisVer[Inidx].end(),
              Vertex.LoopBasis[IN].begin());
    std::copy(LoopBasisVer[Outidx].begin(), LoopBasisVer[Outidx].end(),
              Vertex.LoopBasis[OUT].begin());
    VerIndex.push_back(_AddOneVerToPool(Pool, Vertex));
  }
  return VerIndex;
}

vector<vertex4 *> _AddAllVer4ToPool(pool &Pool, vector<tau> &VerBasis,
                                    vector<int> &Ver4Legs,
                                    vector<loop> &LoopBasisG,
                                    vector<loop> &LoopBasisVer,
                                    vector<int> &VerType, int Ver4Num) {
  vector<vertex4 *> Ver4Index;
  for (int i = 0; i < Ver4Num; i++) {
    int Inidx = 2 * i, Outidx = 2 * i + 1;
    // construct a new vertex 4 function
    vertex4 Vertex4;
    Vertex4.Type = {VerType[Inidx], VerType[Outidx]};

    // build Interaction loop basis
    Vertex4.IntLoopBasis[DIRECT].fill(0);
    Vertex4.IntLoopBasis[EXCHANGE].fill(0);
    std::copy(LoopBasisVer[Inidx].begin(), LoopBasisVer[Inidx].end(),
              Vertex4.IntLoopBasis[DIRECT].begin());
    std::copy(LoopBasisVer[Outidx].begin(), LoopBasisVer[Outidx].end(),
              Vertex4.IntLoopBasis[EXCHANGE].begin());

    // build 4-leg loop basis
    for (int leg = 0; leg < 4; leg++) {
      int legidx = 4 * i + leg; // index shift
      Vertex4.LoopBasis[leg].fill(0);
      int gidx = Ver4Legs[legidx];
      std::copy(LoopBasisG[gidx].begin(), LoopBasisG[gidx].end(),
                Vertex4.LoopBasis[leg].begin());
    }

    if (Para.UseVer4)
      Ver4Index.push_back(_AddOneVer4ToPool(Pool, Vertex4));
    else
      Ver4Index.push_back(_AddOneInteractionToPool(Pool, Vertex4));
  }
  return Ver4Index;
}

diagram ReadOneDiagram(istream &DiagFile, pool &Pool, int Order, int LoopNum,
                       int GNum, int Ver4Num) {
  string buff;
  diagram Diagram;

  Diagram.G.fill(nullptr);
  Diagram.Ver.fill(nullptr);
  Diagram.Ver4.fill(nullptr);

  //////// Diagram Topology  ////////////////////////
  buff = _GetOneLine(DiagFile); // title
  vector<int> Permutation = _ExtractOneLine<int>(DiagFile);
  copy(Permutation.begin(), Permutation.end(), Diagram.Permutation.begin());

  //////// symmetry factor //////////////////
  buff = _GetOneLine(DiagFile); // title
  Diagram.SymFactor = _ExtractOneLine<double>(DiagFile)[0];

  //////// Propagator type //////////////////
  buff = _GetOneLine(DiagFile); // title
  auto GType = _ExtractOneLine<int>(DiagFile);

  /////// Ver Basis  /////////////////////////
  buff = _GetOneLine(DiagFile); // title
  vector<int> StartVer = _ExtractOneLine<int>(DiagFile);
  vector<int> EndVer = _ExtractOneLine<int>(DiagFile);
  vector<tau> VerBasis;
  for (int i = 0; i < GNum; i++)
    VerBasis.push_back(tau({StartVer[i], EndVer[i]}));

  /////// Loop Basis  /////////////////////////
  buff = _GetOneLine(DiagFile); // title
  vector<vector<double>> TransposedLoopBasis;
  for (int j = 0; j < LoopNum; j++)
    TransposedLoopBasis.push_back(_ExtractOneLine<double>(DiagFile));

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
  vector<vector<double>> TransposedLoopBasisVer;
  if (Ver4Num > 0)
    for (int j = 0; j < LoopNum; j++)
      TransposedLoopBasisVer.push_back(_ExtractOneLine<double>(DiagFile));
  vector<loop> LoopBasisVer = _Transpose(TransposedLoopBasisVer);

  /////// Spin Factor  ////////////////////
  buff = _GetOneLine(DiagFile); // title
  auto SpinFactor = _ExtractOneLine<double>(DiagFile);
  copy(SpinFactor.begin(), SpinFactor.end(),
       Diagram.SpinFactor.begin()); // copy spin factor into the group member

  ////////   Add G to GPool ////////////////////////
  vector<green *> GIndex =
      _AddAllGToPool(Pool, VerBasis, LoopBasis, GType, GNum);
  ASSERT_ALLWAYS(GIndex.size() == GNum,
                 "Number of Green's function does not match!");
  copy(GIndex.begin(), GIndex.end(), Diagram.G.begin());

  //   cout << "in diagram" << endl;
  //   cout << ToString(*(Diagram.GIndex[0])) << endl;
  //   cout << Diagram.GIndex[0] << endl;

  if (Ver4Num > 0) {
    //////  Add 4-Vertex to Ver4Pool /////////
    vector<vertex4 *> Ver4Index = _AddAllVer4ToPool(
        Pool, VerBasis, Ver4Legs, LoopBasis, LoopBasisVer, VerType, Ver4Num);
    copy(Ver4Index.begin(), Ver4Index.end(), Diagram.Ver4.begin());
    //////  Add Vertex to VerPool /////////
    vector<vertex *> VerIndex =
        _AddAllVerToPool(Pool, VerBasis, LoopBasisVer, VerType, Ver4Num);
    copy(VerIndex.begin(), VerIndex.end(), Diagram.Ver.begin());
  }

  return Diagram;
}

group diag::ReadOneGroup(istream &DiagFile, pool &Pool) {
  group Group;
  string buff = _GetOneLine(DiagFile); // group type, simply skip

  Group.HugenNum = _ExtractOneLine<int>(DiagFile)[0];
  ASSERT_ALLWAYS(Group.HugenNum <= MaxDiagNum,
                 "Diagram Number must be smaller than " << MaxDiagNum);

  Group.Order = _ExtractOneLine<int>(DiagFile)[0];
  ASSERT_ALLWAYS(Group.Order <= MaxOrder,
                 "Order Number must be smaller than " << MaxOrder);

  Group.GNum = _ExtractOneLine<int>(DiagFile)[0];
  ASSERT_ALLWAYS(Group.GNum <= MaxGNum,
                 "G Number must be smaller than " << MaxGNum);

  Group.Ver4Num = _ExtractOneLine<int>(DiagFile)[0];
  ASSERT_ALLWAYS(Group.Ver4Num <= MaxVer4Num,
                 "Ver4 Number must be smaller than " << MaxVer4Num);

  Group.LoopNum = _ExtractOneLine<int>(DiagFile)[0];
  ASSERT_ALLWAYS(Group.LoopNum <= MaxLoopNum,
                 "Loop Number must be smaller than " << MaxLoopNum);

  vector<int> ExtLoop = _ExtractOneLine<int>(DiagFile);
  Group.ExtLoopNum = ExtLoop.size();
  Group.IsExtLoop.fill(false);
  for (auto index : ExtLoop)
    Group.IsExtLoop[index] = true;

  vector<int> LockedLoop = _ExtractOneLine<int>(DiagFile);
  Group.IsLockedLoop.fill(false);
  for (auto index : LockedLoop)
    Group.IsLockedLoop[index] = true;

  Group.TauNum = _ExtractOneLine<int>(DiagFile)[0];
  ASSERT_ALLWAYS(Group.TauNum <= MaxTauNum,
                 "Tau Number must be smaller than " << MaxTauNum);

  vector<int> ExtTau = _ExtractOneLine<int>(DiagFile);
  Group.ExtTauNum = ExtTau.size();
  Group.IsExtTau.fill(false);
  for (auto index : ExtTau)
    Group.IsExtTau[index] = true;

  vector<int> LockedTau = _ExtractOneLine<int>(DiagFile);
  Group.IsLockedTau.fill(false);
  for (auto index : LockedTau)
    Group.IsLockedTau[index] = true;

  Group.InternalLoopNum = Group.LoopNum - Group.ExtLoopNum;
  Group.InternalTauNum = Group.TauNum - Group.ExtTauNum;

  for (int i = 0; i < Group.HugenNum; i++) {
    diagram Diagram = ReadOneDiagram(DiagFile, Pool, Group.Order, Group.LoopNum,
                                     Group.GNum, Group.Ver4Num);
    Diagram.ID = i;
    Group.Diag.push_back(Diagram);
    // for (int i = 0; i < Group.GNum; i++)
    //   cout << ToString(*(Group.DiagList.back().GIndex[i])) << endl;
  }
  diag::Test(Group);
  //   cout << "Group" << endl;
  //   for (int i = 0; i < Group.GNum; i++)
  //   cout << "in Group" << endl;
  //   cout << ToString(Group) << endl;
  //   cout << ToString(*(Group.DiagList[0].GIndex[0])) << endl;

  return Group;
}

std::string ToString(const diagram &Diag) {
  std::ostringstream oss;
  oss << "DiagID: " << Diag.ID << "\n Weight=" << Diag.Weight
      << "\n SymFactor=" << Diag.SymFactor << endl;
  return oss.str();
}
std::string ToString(const group &Group) {
  return fmt::format("GroupID: {0}\n Weight={1}\n NewWeight={2}\n "
                     "HugenNum={3}\n ReWeight={4}\n",
                     Group.ID, Group.Weight, Group.NewWeight, Group.HugenNum,
                     Group.ReWeight);
};
std::string ToString(const green &G) {
  std::ostringstream oss;
  oss << "Version: " << G.Version << "\n Type=" << G.Type
      << "\n Weight=" << G.Weight
      << "\n TauBasis=" << ToString<int>(G.TauBasis.data(), 2)
      << "\n LoopBasis="
      << ToString<double>(G.LoopBasis.data(), (size_t)G.LoopBasis.size())
      << endl;
  return oss.str();
};
std::string ToString(const vertex &);
std::string ToString(const vertex4 &);

void diag::Test(group &Group) {
  for (auto &d : Group.Diag) {
    for (auto i = 0; i < Group.GNum; i++) {
      CHECKNULL(d.G[i]);
    }
    for (auto i = 0; i < Group.Ver4Num; i++) {
      CHECKNULL(d.Ver[i]);
      CHECKNULL(d.Ver4[i]);
    }
  }
}
