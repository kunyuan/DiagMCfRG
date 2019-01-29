#include "weight.h"
#include "utility/abort.h"
#include "utility/vector.h"
#include <array>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>

using namespace diag;
using namespace std;

void weight::ReadDiagrams() {
  Pool.GPoolSize = 0;
  Pool.VerPoolSize = 0;
  Pool.Ver4PoolSize = 0;

  for (int &id : Para.GroupID) {
    // construct filename based on format string and group id
    string FileName = Format(Para.DiagFileFormat, id);
    ifstream DiagFile(FileName);
    ASSERT_ALLWAYS(DiagFile.is_open(),
                   "Unable to find the file " << FileName << endl);
    // group Group;
    LOG_INFO("Find " << FileName << "\n");
    // vector<green> GList;
    istream &DiagFileStream = DiagFile;
    Groups.push_back(ReadOneGroup(DiagFileStream, Pool));
    Groups.back().ID = id;
  }
  LOG_INFO("Find " << Pool.GPoolSize << " indepdent green's function.");
  LOG_INFO("Find " << Pool.VerPoolSize << " indepdent interactions.");
  LOG_INFO("Find " << Pool.Ver4PoolSize << " indepdent 4-vertex.");

  // cout << "After read" << endl;
  // cout << ToString(*(GroupList[0].DiagList[0].GIndex[0])) << endl;
  Initialization();
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
      mom[i] = Random.urn() * Para.Kf / sqrt(D);

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
    // the external momentum only has x component
    Var.ExtMomTable[i][0] = i * Para.MaxExtMom / ExtMomBinSize;
    for (int j = 1; j < D; j++)
      Var.ExtMomTable[i][j] = 0.0;
  }
  Var.CurrExtMomBin = 0;
  // Var.LoopMom[0].fill(0.0);
  // for (int i = 0; i < D; i++)
  //   Var.LoopMom[0][i] = Var.ExtMomTable[Var.CurrExtMomBin][i];
  Var.LoopMom[0] = Var.ExtMomTable[Var.CurrExtMomBin];

  // initialize external tau
  Var.Tau[0] = 0.0;
  Var.Tau[1] = 1.0e-10; // do not make Tau[1]==Tau[0], otherwise the Green's
                        // function is not well-defined
  Var.CurrTau = Var.Tau[1] - Var.Tau[0];

  // initialize group

  Var.CurrVersion = 0;
  //   Var.CurrGroup = &Groups[0];

  Var.CurrGroup = &Groups[0];
  LOG_INFO("Calculating the weights of all objects...")

  ChangeGroup(*Var.CurrGroup, true);
  GetNewWeight(*Var.CurrGroup);
  AcceptChange(*Var.CurrGroup);

  LOG_INFO("Initializating variables done.")
}

void weight::ChangeGroup(group &Group, bool Forced) {
  // the objects (G, Ver or Ver4) in the new group will be recalculated if the
  // either of the following conditions is met: 1) Forced=true, then all objects
  // are forced recalculated 2) object.Version<CurrVersion, means the objects
  // are not in the current group, and are already outdated
  for (auto &d : Group.Diag) {
    // cout << "diag ID: " << d.ID << endl;
    for (int i = 0; i < Group.GNum; i++) {
      // cout << "G: " << i << endl;
      green *G = d.G[i];
      if (Forced || G->Version < Var.CurrVersion) {
        double Tau = Var.Tau[G->TauBasis[OUT]] - Var.Tau[G->TauBasis[IN]];
        G->Excited = true;
        GetMom(G->LoopBasis, Group.LoopNum);
        G->NewWeight = Green(Tau, _Mom, UP, G->Type);
      }
    }
    for (int i = 0; i < Group.Ver4Num; i++) {
      // cout << "Ver: " << i << endl;
      if (UseVertex4) {
        THROW_ERROR(NOTIMPLEMENTED, "Vertex4 not implemented!");
        // vertex4 *Ver4 = d.Ver4Index[i];
        // Ver4->Excited=true;
        // Ver4->NewWeight=
      } else {
        vertex *Ver = d.Ver[i];
        if (Forced || Ver->Version < Var.CurrVersion) {
          Ver->Excited = {true, true};
          if (!IsInteractionReducible(Ver->LoopBasis[IN], Group.LoopNum)) {
            GetMom(Ver->LoopBasis[IN], Group.LoopNum);
            Ver->NewWeight[IN] = Interaction(0.0, _Mom, Ver->Type[IN]);
          } else {
            Ver->NewWeight[IN] = 0.0;
          }
          if (!IsInteractionReducible(Ver->LoopBasis[OUT], Group.LoopNum)) {
            GetMom(Ver->LoopBasis[OUT], Group.LoopNum);
            Ver->NewWeight[OUT] = Interaction(0.0, _Mom, Ver->Type[OUT]);
          } else {
            Ver->NewWeight[OUT] = 0.0;
          }
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
        GetMom(G->LoopBasis, Group.LoopNum);
        G->NewWeight = Green(Tau, _Mom, UP, G->Type);
      }
    }
    for (int i = 0; i < Group.Ver4Num; i++) {
      if (UseVertex4) {
        THROW_ERROR(NOTIMPLEMENTED, "Vertex4 not implemented!");
      } else {
        vertex *Ver = d.Ver[i];
        if (Ver->LoopBasis[IN][MomIndex] != 0) {
          Ver->Excited[IN] = true;
          if (!IsInteractionReducible(Ver->LoopBasis[IN], Group.LoopNum)) {
            GetMom(Ver->LoopBasis[IN], Group.LoopNum);
            Ver->NewWeight[IN] = Interaction(0.0, _Mom, Ver->Type[IN]);
          } else {
            Ver->NewWeight[IN] = 0.0;
          }
        }
        if (Ver->LoopBasis[OUT][MomIndex] != 0) {
          Ver->Excited[OUT] = true;
          if (!IsInteractionReducible(Ver->LoopBasis[OUT], Group.LoopNum)) {
            GetMom(Ver->LoopBasis[OUT], Group.LoopNum);
            Ver->NewWeight[OUT] = Interaction(0.0, _Mom, Ver->Type[OUT]);
          } else {
            Ver->NewWeight[OUT] = 0.0;
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
      // we also assume that 0 will never be changed!
      if (TauIndex <= 1 && (TauIn == 1 || TauOut == 1))
        ReCalcFlag = true;
      // if TauIndex is >1, then G is recaculated if TauIndex/2==TauIn/2 or
      // TauOut/2
      if (TauIndex > 1 &&
          (TauIn / 2 == TauIndex / 2 || TauOut / 2 == TauIndex / 2))
        ReCalcFlag = true;

      if (ReCalcFlag) {
        // trigger recalculation
        double Tau = Var.Tau[TauOut] - Var.Tau[TauIn];
        G->Excited = true;
        GetMom(G->LoopBasis, Group.LoopNum);
        G->NewWeight = Green(Tau, _Mom, UP, G->Type);
      }
    }
  }
}

double weight::GetNewWeight(group &Group) {
  static double VIn, VOut, TotWeight;
  Group.NewWeight = 0.0;

  for (auto &d : Group.Diag) {
    double GWeight = d.SymFactor;
    for (int i = 0; i < Group.GNum; i++) {
      if (d.G[i]->Excited)
        GWeight *= d.G[i]->NewWeight;
      else
        GWeight *= d.G[i]->Weight;
    }

    double VerWeight;
    if (UseVertex4) {
      THROW_ERROR(NOTIMPLEMENTED, "Ver4 has not yet been implemented!");
    } else {
      if (Group.Ver4Num == 0) {
        VerWeight = d.SpinFactor[0];
        // cout << "spin factor: " << d.SpinFactor[0] << endl;
      } else {

        vertex *Ver = d.Ver[0];
        _Tree[0][0] = Ver->Excited[IN] ? Ver->NewWeight[IN] : Ver->Weight[IN];
        _Tree[0][1] =
            Ver->Excited[OUT] ? Ver->NewWeight[OUT] : Ver->Weight[OUT];

        int BlockNum = 2;
        for (int level = 1; level < Group.Ver4Num; level++) {

          vertex *Ver = d.Ver[level];
          VIn = Ver->Excited[IN] ? Ver->NewWeight[IN] : Ver->Weight[IN];
          VOut = Ver->Excited[OUT] ? Ver->NewWeight[OUT] : Ver->Weight[OUT];

          for (int j = 0; j < BlockNum; j++) {
            _Tree[level][2 * j] = _Tree[level - 1][j] * VIn;
            _Tree[level][2 * j + 1] = _Tree[level - 1][j] * VOut;
          }
          BlockNum *= 2;
        }

        // cout << BlockNum << endl;
        // cout << d.SpinFactor[0] << ", " << d.SpinFactor[1] << endl;

        VerWeight = 0.0;
        for (int j = 0; j < BlockNum; j++)
          VerWeight += _Tree[Group.Ver4Num - 1][j] * d.SpinFactor[j];

        //============= for spin case ===========================//

        // double TempWeightIn, TempWeightOut;
        // for (int i = 0; i < Group.Ver4Num; i++) {
        //   //=========== for spinless case ===========================//
        //   vertex *Ver = d.Ver[i];
        //   TempWeightIn =
        //       Ver->Excited[IN] ? Ver->NewWeight[IN] : Ver->Weight[IN];
        //   TempWeightOut =
        //       Ver->Excited[OUT] ? Ver->NewWeight[OUT] : Ver->Weight[OUT];
        //   VerWeight *= TempWeightIn - TempWeightOut;
      }
    }
    d.NewWeight = GWeight * VerWeight / pow(2 * PI, D * Group.InternalLoopNum);

    // Group Weight= sum of all diagram weight in the group
    Group.NewWeight += d.NewWeight;
  }
  return Group.NewWeight;
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

void weight::RejectChange(group &Group) {
  for (auto &d : Group.Diag) {
    for (int i = 0; i < Group.GNum; i++) {
      if (d.G[i]->Excited)
        d.G[i]->Excited = false;
      for (int i = 0; i < Group.Ver4Num; i++) {
        if (UseVertex4) {
          THROW_ERROR(NOTIMPLEMENTED, "Ver4 has not yet been implemented!");
        } else {
          if (d.Ver[i]->Excited[0])
            d.Ver[i]->Excited[0] = false;
          if (d.Ver[i]->Excited[1])
            d.Ver[i]->Excited[1] = false;
        }
      }
    }
  }
}

void weight::GetMom(const loop &LoopBasis, const int &LoopNum) {
  // In C++11, because of the move semantics, there is no additional cost by
  // returning an array

  auto &loopmom = Var.LoopMom;
  for (int d = 0; d < D; ++d)
    _Mom[d] = loopmom[0][d] * LoopBasis[0];

  for (int i = 1; i < LoopNum; ++i)
    for (int d = 0; d < D; ++d)
      _Mom[d] += loopmom[i][d] * LoopBasis[i];
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