#include "weight.h"
#include "utility/abort.h"
#include "utility/vector.h"
#include <array>
#include <iostream>
#include <sstream>
#include <string>

using namespace diag;
using namespace vertex;
using namespace std;

void weight::ReadDiagrams(string FilePrefix) {
  Pool.GPoolSize = 0;
  Pool.VerPoolSize = 0;
  Pool.Ver4PoolSize = 0;

  int GroupID = 0;
  while (true) {
    GroupID++;
    ifstream DiagFile(FilePrefix + to_string(GroupID) + ".txt");
    if (DiagFile) {
      // group Group;
      LOG_INFO("Find " << FilePrefix + to_string(GroupID) + ".txt\n");
      // vector<green> GList;
      istream &DiagFileStream = DiagFile;
      Groups.push_back(ReadOneGroup(DiagFileStream, Pool));
      Groups.back().ID = GroupID;
      // ReadOneGroup(DiagFile, Group, GList);
      // cout << "OutGroup" << endl;
      // cout << ToString(*(Group.DiagList[0].GIndex[0])) << endl;
      // cout << ToString(*(GroupList[0].DiagList[0].GIndex[0])) << endl;
      // cout << ToString(Group) << endl;
      // cout << ToString(GroupList[0]) << endl;
      // cout << ToString(Pool.GPool[0]) << endl;
    } else
      break;
  }
  LOG_INFO("Find " << Pool.GPoolSize << " indepdent green's function.");
  LOG_INFO("Find " << Pool.VerPoolSize << " indepdent interactions.");
  LOG_INFO("Find " << Pool.Ver4PoolSize << " indepdent 4-vertex.");

  // cout << "After read" << endl;
  // cout << ToString(*(GroupList[0].DiagList[0].GIndex[0])) << endl;
}

void weight::Initialization() {

  LOG_INFO("Initializating diagram states ...")
  for (auto &group : Groups) {
    group.ReWeight = 1.0;
    for (auto &diag : group.Diag) {
      for (int i = 0; i < group.GNum; i++) {
        diag.G[i]->Excited = false;
        diag.G[i]->Version = -1;
        diag.G[i]->Weight = 1.0e-10;
      }
      for (int i = 0; i < group.Ver4Num; i++) {
        diag.Ver[i]->Excited = {false, false};
        diag.Ver[i]->Version = -1;
        diag.Ver[i]->Weight = {1.0e-10, 1.0e-10};
      }
      for (int i = 0; i < group.Ver4Num; i++) {
        diag.Ver4[i]->Excited = false;
        diag.Ver4[i]->Version = -1;
        diag.Ver4[i]->Weight = 1.0e-10;
      }
    }
  }

  LOG_INFO("Initializating MC variables ...")
  // initialize momentum variables
  for (auto &mom : Var.LoopMom)
    for (int i = 0; i < D; i++)
      mom[i] = Random.urn() * Para.Kf;

  // initialize tau variables
  for (int i = 0; i < MaxTauNum / 2; i++) {
    Var.Tau[2 * i] = Random.urn() * Para.Beta;
    Var.Tau[2 * i + 1] = Var.Tau[2 * i]; // assume even and odd tau are the same
  }

  // initialize spin variables
  for (auto &sp : Var.LoopSpin)
    sp = (spin)(Random.irn(0, 1));

  // initialize external momentum
  for (int i = 0; i < ExtMomBinSize; i++) {
    Var.ExtMomTable[i] = i * Para.MaxExtMom / ExtMomBinSize;
  }
  Var.CurrExtMomBin = 0;
  Var.LoopMom[0].fill(0.0);
  Var.LoopMom[0][0] = Var.ExtMomTable[Var.CurrExtMomBin];

  // initialize external tau
  Var.Tau[0] = 0.0;
  Var.CurrTau = Var.Tau[1] - Var.Tau[0];

  // initialize group

  Var.CurrVersion = 0;
  Var.CurrGroup = &Groups[0];
  LOG_INFO("Calculating the weights of all objects...")

  ChangeGroup(*Var.CurrGroup, true);
  Var.CurrWeight = GetNewWeight(*Var.CurrGroup);
  AcceptChange(*Var.CurrGroup);

  LOG_INFO("Initializating variables done.")
}

double weight::GetNewWeight(group &Group) {
  Group.NewWeight = 1.0;
  for (auto &d : Group.Diag) {
    double GWeight = d.SymFactor;
    for (int i = 0; i < Group.GNum; i++) {
      if (d.G[i]->Excited)
        GWeight *= d.G[i]->NewWeight;
      else
        GWeight *= d.G[i]->Weight;
    }
    double VerWeight = 1.0;
    double TempWeightIn, TempWeightOut;
    for (int i = 0; i < Group.Ver4Num; i++)
      if (UseVertex4) {
        THROW_ERROR(NOTIMPLEMENTED, "Ver4 has not yet been implemented!");
      } else {
        // only for spinless case
        vertex *Ver = d.Ver[i];
        TempWeightIn = Ver->Excited[IN] ? Ver->NewWeight[IN] : Ver->Weight[IN];
        TempWeightOut =
            Ver->Excited[OUT] ? Ver->NewWeight[OUT] : Ver->Weight[OUT];
        VerWeight *= TempWeightIn - TempWeightOut;
      }
    d.NewWeight =
        GWeight * VerWeight / pow(1.0 / 2 / PI, D * Group.InternalLoopNum);
    Group.NewWeight *= d.NewWeight;
  }
  return Group.NewWeight * Group.ReWeight;
}

void weight::AcceptChange(group &Group) {
  Var.CurrVersion++;
  Var.CurrGroup = &Group;
  Group.Weight = Group.NewWeight; // accept group  newweight
  for (auto &d : Group.Diag) {
    d.Weight = d.NewWeight; // accept diagram newweight
    for (int i = 0; i < Group.GNum; i++) {
      green *G = d.G[i];
      G->Version = Var.CurrVersion;
      if (G->Excited) {
        G->Excited = false;
        G->Weight = G->NewWeight;
      }
    }
    for (int i = 0; i < Group.Ver4Num; i++)
      if (UseVertex4) {
        THROW_ERROR(NOTIMPLEMENTED, "Ver4 has not yet been implemented!");
      } else {
        vertex *Ver = d.Ver[i];
        Ver->Version = Var.CurrVersion;
        if (Ver->Excited[IN]) {
          Ver->Excited[IN] = false;
          Ver->Weight[IN] = Ver->NewWeight[IN];
        }
        if (Ver->Excited[OUT]) {
          Ver->Excited[OUT] = false;
          Ver->Weight[OUT] = Ver->NewWeight[OUT];
        }
      }
  }
}

momentum weight::_Mom(const loop &LoopBasis, const int &LoopNum) {
  // In C++11, because of the move semantics, there is no additional cost by
  // returning an array
  momentum Mom;
  for (int d = 0; d < D; d++)
    for (int i = 0; i < LoopNum; i++)
      Mom[d] += Var.LoopMom[i][d] * LoopBasis[i];
  return Mom;
}

void weight::ChangeGroup(group &Group, bool Forced) {
  // the objects (G, Ver or Ver4) in the new group will be recalculated if the
  // either of the following conditions is met: 1) Forced=true, then all objects
  // are forced recalculated 2) object.Version<CurrVersion, means the objects
  // are not in the current group, and are already outdated
  for (auto &d : Group.Diag) {
    cout << "diag ID: " << d.ID << endl;
    for (int i = 0; i < Group.GNum; i++) {
      cout << "G: " << i << endl;
      green *G = d.G[i];
      if (Forced || G->Version < Var.CurrVersion) {
        double Tau = Var.Tau[G->TauBasis[OUT]] - Var.Tau[G->TauBasis[IN]];
        G->Excited = true;
        G->NewWeight =
            Green(Tau, _Mom(G->LoopBasis, Group.LoopNum), UP, G->Type);
      }
    }
    for (int i = 0; i < Group.Ver4Num; i++) {
      cout << "Ver: " << i << endl;
      if (UseVertex4) {
        THROW_ERROR(NOTIMPLEMENTED, "Vertex4 not implemented!");
        // vertex4 *Ver4 = d.Ver4Index[i];
        // Ver4->Excited=true;
        // Ver4->NewWeight=
      } else {
        vertex *Ver = d.Ver[i];
        if (Forced || Ver->Version < Var.CurrVersion) {
          Ver->Excited = {true, true};
          Ver->NewWeight[0] = Interaction(
              0.0, _Mom(Ver->LoopBasis[0], Group.LoopNum), Ver->Type[0]);
          Ver->NewWeight[1] = Interaction(
              0.0, _Mom(Ver->LoopBasis[1], Group.LoopNum), Ver->Type[1]);
        }
      }
    }
  }
}

void weight::ChangeMom(group &Group, int MomIndex) {
  for (auto &d : Group.Diag) {
    for (int i = 0; i < Group.GNum; i++) {
      green *G = d.G[i];
      if (G->LoopBasis[MomIndex] != 0) {
        double Tau = Var.Tau[G->TauBasis[OUT]] - Var.Tau[G->TauBasis[IN]];
        G->Excited = true;
        G->NewWeight =
            Green(Tau, _Mom(G->LoopBasis, Group.LoopNum), UP, G->Type);
      }
    }
    for (int i = 0; i < Group.Ver4Num; i++) {
      if (UseVertex4) {
        THROW_ERROR(NOTIMPLEMENTED, "Vertex4 not implemented!");
      } else {
        vertex *Ver = d.Ver[i];
        if (Ver->LoopBasis[0][MomIndex] != 0) {
          Ver->Excited[0] = true;
          if (IsInteractionReducible(Ver->LoopBasis[0], Group.LoopNum)) {
            Ver->NewWeight[0] = Interaction(
                0.0, _Mom(Ver->LoopBasis[0], Group.LoopNum), Ver->Type[0]);
          } else {
            Ver->NewWeight[0] = 0.0;
          }
        }
        if (Ver->LoopBasis[1][MomIndex] != 0) {
          Ver->Excited[1] = true;
          if (IsInteractionReducible(Ver->LoopBasis[0], Group.LoopNum)) {
            Ver->NewWeight[1] = Interaction(
                0.0, _Mom(Ver->LoopBasis[1], Group.LoopNum), Ver->Type[1]);
          } else {
            Ver->NewWeight[1] = 0.0;
          }
        }
      }
    }
  }
}

void weight::ChangeTau(group &Group, int TauIndex) {
  // TODO: we assume TauLeft==TauRight for now
  for (auto &d : Group.Diag) {
    for (int i = 0; i < Group.GNum; i++) {
      green *G = d.G[i];
      int TauIn = G->TauBasis[IN];
      int TauOut = G->TauBasis[OUT];
      bool ReCalcFlag = false;

      // if TauIndex is 0, or 1, then G is recaculated only if TauIn or TauOut
      // exactly matches TauIndex
      if (TauIndex <= 1 && (TauIn == TauIndex || TauOut == TauIndex))
        ReCalcFlag = true;
      // if TauIndex is >1, then G is recaculated if TauIndex/2==TauIn/2 or
      // TauOut/2
      if (TauIndex > 1 && (TauIn / 2 == TauIndex / 2 || TauOut == TauIndex / 2))
        ReCalcFlag = true;

      if (ReCalcFlag) {
        // trigger recalculation
        double Tau = Var.Tau[TauOut] - Var.Tau[TauIn];
        G->Excited = true;
        G->NewWeight =
            Green(Tau, _Mom(G->LoopBasis, Group.LoopNum), UP, G->Type);
      }
    }
  }
}

bool weight::IsInteractionReducible(loop &LoopBasisVer, int LoopNum) {
  // check if an interaction is reducible
  if (abs(LoopBasisVer[0]) != 1)
    return false;
  bool Flag = true;
  for (int i = 1; i < LoopNum; i++) {
    if (LoopBasisVer[i] != 0) {
      Flag = false;
      break;
    }
  }
  return Flag;
};