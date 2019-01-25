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
  Pool.GPool.clear();
  Pool.VerPool.clear();
  Pool.Ver4Pool.clear();

  int GroupID = 0;
  while (true) {
    GroupID++;
    ifstream DiagFile(FilePrefix + to_string(GroupID) + ".txt");
    if (DiagFile) {
      // group Group;
      LOG_INFO("Find " << FilePrefix + to_string(GroupID) + ".txt\n");
      // vector<green> GList;
      istream &DiagFileStream = DiagFile;
      group Group = ReadOneGroup(DiagFileStream, Pool);
      Group.ID = GroupID;
      GroupList.push_back(Group);
      // ReadOneGroup(DiagFile, Group, GList);
      cout << "OutGroup" << endl;
      cout << ToString(*(Group.DiagList[0].GIndex[0])) << endl;
      cout << ToString(*(GroupList[0].DiagList[0].GIndex[0])) << endl;
      cout << ToString(Group) << endl;
      cout << ToString(GroupList[0]) << endl;
      cout << ToString(Pool.GPool[0]) << endl;
    } else
      break;
  }
  LOG_INFO("Find " << Pool.GPool.size() << " indepdent green's function.");
  LOG_INFO("Find " << Pool.VerPool.size() << " indepdent interactions.");
  LOG_INFO("Find " << Pool.Ver4Pool.size() << " indepdent 4-vertex.");

  cout << "After read" << endl;
  cout << ToString(*(GroupList[0].DiagList[0].GIndex[0])) << endl;
}

void weight::Initialization() {
  LOG_INFO("Initializating variables ...")

  Variable.CurrVersion = 0;
  int GroupIndex = 0;
  Variable.CurrGroup = &(GroupList[GroupIndex]);
  ChangeGroup(GroupList[GroupIndex], true);

  for (auto &mom : Variable.LoopMom)
    for (int i = 0; i < D; i++)
      mom[i] = Random.urn();

  for (int i = 0; i < MaxTauNum / 2; i++) {
    Variable.Tau[2 * i] = Random.urn();
    Variable.Tau[2 * i + 1] = Random.urn();
  }

  for (auto &sp : Variable.LoopSpin)
    sp = (spin)(Random.irn(0, 1));

  LOG_INFO("Calculating the weights of all objects...")
  Variable.CurrWeight = GetNewWeight(*Variable.CurrGroup);
  LOG_INFO("Initializating variables done.")
}

double weight::GetNewWeight(group &Group) {
  Group.NewWeight = 1.0;
  for (auto &d : Group.DiagList) {
    double GWeight = d.SymFactor;
    for (int i = 0; i < Group.GNum; i++) {
      if (d.GIndex[i]->Excited)
        GWeight *= d.GIndex[i]->NewWeight;
      else
        GWeight *= d.GIndex[i]->Weight;
    }
    double VerWeight = 1.0;
    double TempWeightIn, TempWeightOut;
    for (int i = 0; i < Group.Ver4Num; i++)
      if (UseVertex4) {
        THROW_ERROR(NOTIMPLEMENTED, "Ver4 has not yet been implemented!");
      } else {
        // only for spinless case
        vertex *Ver = d.VerIndex[i];
        TempWeightIn = Ver->Excited[IN] ? Ver->NewWeight[IN] : Ver->Weight[IN];
        TempWeightOut =
            Ver->Excited[OUT] ? Ver->NewWeight[OUT] : Ver->Weight[OUT];
        VerWeight *= TempWeightIn - TempWeightOut;
      }
    d.NewWeight =
        GWeight * VerWeight / pow(1.0 / 2 / PI, D * Group.InternalLoopNum);
    Group.NewWeight *= d.NewWeight;
  }
  return Group.NewWeight;
}

void weight::AcceptChange(group &Group) {
  Variable.CurrVersion++;
  Variable.CurrGroup = &Group;
  Group.Weight = Group.NewWeight;
  for (auto &d : Group.DiagList) {
    d.Weight = d.NewWeight;
    for (int i = 0; i < Group.GNum; i++) {
      green *G = d.GIndex[i];
      G->Version = Variable.CurrVersion;
      if (G->Excited) {
        G->Excited = false;
        G->Weight = G->NewWeight;
      }
    }
    for (int i = 0; i < Group.Ver4Num; i++)
      if (UseVertex4) {
        THROW_ERROR(NOTIMPLEMENTED, "Ver4 has not yet been implemented!");
      } else {
        vertex *Ver = d.VerIndex[i];
        Ver->Version = Variable.CurrVersion;
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
      Mom[d] += Variable.LoopMom[i][d] * LoopBasis[i];
  return Mom;
}

void weight::ChangeGroup(group &Group, bool Forced) {
  // the objects (G, Ver or Ver4) in the new group will be recalculated if the
  // either of the following conditions is met: 1) Forced=true, then all objects
  // are forced recalculated 2) object.Version<CurrVersion, means the objects
  // are not in the current group, and are already outdated
  for (auto &d : Group.DiagList) {
    cout << "diag ID: " << d.ID << endl;
    for (int i = 0; i < Group.GNum; i++) {
      cout << "G: " << i << endl;
      green *G = d.GIndex[i];
      if (Forced || G->Version < Variable.CurrVersion) {
        double Tau =
            Variable.Tau[G->TauBasis[OUT]] - Variable.Tau[G->TauBasis[IN]];
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
        vertex *Ver = d.VerIndex[i];
        if (Forced || Ver->Version < Variable.CurrVersion) {
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

void weight::ChangeTau(group &Group, int TauIndex) {
  // TODO: we assume TauLeft==TauRight for now
  for (auto &d : Group.DiagList) {
    for (int i = 0; i < Group.GNum; i++) {
      green *G = d.GIndex[i];
      int TauIn = G->TauBasis[IN];
      int TauOut = G->TauBasis[OUT];
      if (TauIn / 2 == TauIndex / 2 || TauOut == TauIndex / 2) {
        // if TauIndex is equal to either TauIn or TauOut, trigger recalculation
        double Tau = Variable.Tau[TauOut] - Variable.Tau[TauIn];
        G->Excited = true;
        G->NewWeight =
            Green(Tau, _Mom(G->LoopBasis, Group.LoopNum), UP, G->Type);
      }
    }
  }
}