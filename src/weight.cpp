#include "weight.h"
#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/vector.h"
#include <array>
#include <iostream>
#include <string>

using namespace diag;
using namespace std;

void weight::ReadDiagrams() {
  Pool.GPoolSize = 0;
  Pool.Ver4PoolSize = 0;

  int ID = 0;
  for (auto &name : Para.GroupName) {
    // construct filename based on format string and group id
    string FileName = fmt::format(Para.DiagFileFormat, name);
    ifstream DiagFile(FileName);
    ASSERT_ALLWAYS(DiagFile.is_open(),
                   "Unable to find the file " << FileName << endl);
    // group Group;
    LOG_INFO("Find " << FileName << "\n");
    // vector<green> GList;
    istream &DiagFileStream = DiagFile;
    Groups.push_back(ReadOneGroup(DiagFileStream, Pool));
    Groups.back().Name = name;
    Groups.back().ID = ID;
    ID++;
  }
  LOG_INFO("Find " << Pool.GPoolSize << " indepdent green's function.");
  LOG_INFO("Find " << Pool.Ver4PoolSize << " indepdent 4-vertex.");

  // cout << "After read" << endl;
  // cout << ToString(*(GroupList[0].DiagList[0].GIndex[0])) << endl;
  Initialization();
}

void weight::Initialization() {

  LOG_INFO("Initializating diagram states ...")
  for (auto &group : Groups) {
    group.ReWeight = 1.0;
    // for (auto &diag : group.Diag) {
    //   for (int i = 0; i < group.GNum; i++) {
    //     diag.G[i]->Excited = false;
    //     diag.G[i]->Version = -1;
    //     diag.G[i]->Weight = 1.0e-10;
    //   }
    //   for (int i = 0; i < group.Ver4Num; i++) {
    //     diag.Ver4[i]->Excited = {false, false};
    //     diag.Ver4[i]->Version = -1;
    //     diag.Ver4[i]->Weight = {1.0e-10, -1.0e-10};
    //   }
    // }

    // if (Para.ObsType == EQUALTIME) {
    //   for (int i = 0; i < group.Ver4Num; ++i)
    //     if (group.IsExtTau[i])
    //       // to measure equal-time observable, lock all external tau
    //       group.IsLockedTau[i] = true;
    // }
  }

  LOG_INFO("Initializating MC variables ...")
  // initialize momentum variables
  for (auto &mom : Var.LoopMom) {
    for (int i = 0; i < D; i++)
      mom[i] = Random.urn() * Para.Kf / sqrt(D);
    mom[D] = Para.DeltaW * 0.5;
  }

  // initialize spin variables
  for (auto &sp : Var.LoopSpin)
    sp = (spin)(Random.irn(0, 1));

  Var.CurrExtMomBin = 0;
  // Var.LoopMom[0].fill(0.0);
  // for (int i = 0; i < D; i++)
  //   Var.LoopMom[0][i] = Var.ExtMomTable[Var.CurrExtMomBin][i];
  Var.LoopMom[0] = Para.ExtMomTable[Var.CurrExtMomBin];

  // initialize external tau
  // Var.Tau[0] = 0.0;
  // Var.Tau[1] = 1.0e-10; // do not make Tau[1]==Tau[0], otherwise the Green's
  // function is not well-defined

  // initialize group

  Var.CurrVersion = 0;
  //   Var.CurrGroup = &Groups[0];

  Var.CurrGroup = &Groups[0];

  Var.CurrIRScaleBin = ScaleBinSize / 1.5;

  // initialize RG staff
  // Var.CurrScale = ScaleBinSize - 1;
  Var.CurrScale = Para.Kf;

  LOG_INFO("Calculating the weights of all objects...")

  // ChangeGroup(*Var.CurrGroup, true);
  GetNewWeight(*Var.CurrGroup);
  AcceptChange(*Var.CurrGroup);

  LOG_INFO("Initializating variables done.")
}

cmplx weight::GetNewWeight(group &Group) {
  Group.NewWeight = fRG(Group.Order);
  // cout << Group.NewWeight << "order: " << Group.Order << endl;
  return Group.NewWeight;
}

void weight::AcceptChange(group &Group) {
  Var.CurrVersion++;
  Var.CurrGroup = &Group;
  Group.Weight = Group.NewWeight; // accept group  newweight
}

void weight::RejectChange(group &Group) { return; }

void weight::Measure(double WeightFactor) {
  if (Para.Type == RG && Para.Vertex4Type == MOM_ANGLE) {
    if (Var.CurrScale >= Para.ScaleTable[Var.CurrIRScaleBin])
      VerQTheta.Measure(Var.LoopMom[1], Var.LoopMom[2], Var.CurrExtMomBin,
                        Var.CurrScale, Var.CurrGroup->Order, WeightFactor);
  }
}

void weight::Update(double Ratio) {
  if (Para.Type == RG && Para.Vertex4Type == MOM_ANGLE) {
    VerQTheta.Update(Ratio);
  }
}

void weight::Save() {
  if (Para.Type == RG && Para.Vertex4Type == MOM_ANGLE) {
    VerQTheta.Save();
  }
}

void weight::ClearStatis() {
  if (Para.Type == RG && Para.Vertex4Type == MOM_ANGLE) {
    VerQTheta.ClearStatis();
  }
}

void weight::ResetIRScale(int IRScaleBin) {
  Var.CurrIRScaleBin = IRScaleBin;
  if (Para.Type == RG && Para.Vertex4Type == MOM_ANGLE) {
    VerQTheta.ResetIRScale(IRScaleBin);
  }
}

void weight::GetMom(const loop &LoopBasis, const int &LoopNum, momentum &Mom) {
  // In C++11, because of the move semantics, there is no additional cost by
  // returning an array

  auto &loopmom = Var.LoopMom;
  for (int d = 0; d < D; ++d)
    Mom[d] = loopmom[0][d] * LoopBasis[0];

  for (int i = 1; i < LoopNum; ++i)
    for (int d = 0; d < D; ++d)
      Mom[d] += loopmom[i][d] * LoopBasis[i];
}

bool weight::IsInteractionReducible(loop &LoopBasisVer, int LoopNum) {
  // check if an interaction is reducible
  if ((!Equal(LoopBasisVer[0], 1.0)) && (!Equal(LoopBasisVer[0], -1.0)))
    return false;

  bool Flag = true;
  for (int i = 1; i < LoopNum; i++) {
    if (!Equal(LoopBasisVer[i], 0.0)) {
      Flag = false;
      break;
    }
  }
  return Flag;
};

bool weight::IsInteractionReducible(loop &LoopBasisG1, loop &LoopBasisG2,
                                    int LoopNum) {
  // check if an interaction is reducible
  if ((!Equal(LoopBasisG1[0] - LoopBasisG2[0], 1.0)) &&
      (!Equal(LoopBasisG1[0] - LoopBasisG2[0], -1.0)))
    return false;

  bool Flag = true;
  for (int i = 1; i < LoopNum; i++) {
    if (!Equal(LoopBasisG1[i] - LoopBasisG2[i], 0.0)) {
      Flag = false;
      break;
    }
  }
  return Flag;
};
