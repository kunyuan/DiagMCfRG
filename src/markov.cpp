//
//  markov.cpp
//
//  Created by Kun Chen on 1/21/19.
//  Copyright (c) 2019 Kun Chen. All rights reserved.
//
#include "markov.h"

extern parameter Para;
extern RandomFactory Random;

using namespace mc;
using namespace diag;
using namespace std;

#define NAME(x) #x

int markov::DynamicTest() { return Weight.DynamicTest(); }

void markov::Initialization(string FilePrefix) {
  ///==== initialize Weight ============================//
  Weight.ReadDiagrams(FilePrefix);
  Weight.Initialization();

  //===== initialize updates related variable ==========//
  Para.Counter = 0;

  UpdateName[INCREASE_ORDER] = NAME(INCREASE_ORDER);
  UpdateName[DECREASE_ORDER] = NAME(DECREASE_ORDER);
  UpdateName[CHANGE_GROUP] = NAME(CHANGE_GROUP);
  UpdateName[CHANGE_MOM] = NAME(CHANGE_MOM);
  UpdateName[CHANGE_TAU] = NAME(CHANGE_TAU);

  InitialArray(&Accepted[0][0], 1.0e-10, MCUpdates * MaxGroupNum);
  InitialArray(&Proposed[0][0], 1.0e-10, MCUpdates * MaxGroupNum);

  ///==== initialize observable =======================//
  for (auto &g : Weight.Groups) {
    polar p = polar();
    p.fill(1.0e-10);
    Polar.push_back(p);
    PolarStatic.push_back(1.0e-10);
  }
}

void markov::Hop(int Sweep) {
  Para.Counter++;
  for (int i = 0; i < Sweep; i++) {
    double x = Random.urn();
    if (x < 1.0 / MCUpdates)
      IncreaseOrder();
    else if (x < 2.0 / MCUpdates)
      DecreaseOrder();
    else if (x < 3.0 / MCUpdates)
      ChangeGroup();
    else if (x < 4.0 / MCUpdates)
      ChangeMomentum();
    else if (x < 5.0 / MCUpdates)
      ChangeTau();
    else {
    }
  }
}

void markov::PrintMCInfo(){};
void markov::AdjustGroupReWeight(){};

void markov::Measure() {
  double AbsWeight = fabs(Var.CurrWeight);
  double WeightFactor = Var.CurrWeight / AbsWeight * Var.CurrGroup->ReWeight;
  Polar[Var.CurrGroup->ID][Var.CurrExtMomBin] += WeightFactor;
  PolarStatic[Var.CurrGroup->ID] += WeightFactor;
};

void markov::SaveToFile(std::string FilePrefix){};

void markov::IncreaseOrder() {
  group &NewGroup = Groups[Random.irn(0, Groups.size() - 1)];
  if (NewGroup.Order != Var.CurrGroup->Order + 1)
    return;
  Proposed[INCREASE_ORDER][NewGroup.ID] += 1;
};

void markov::DecreaseOrder() {
  group &NewGroup = Groups[Random.irn(0, Groups.size() - 1)];
  if (NewGroup.Order != Var.CurrGroup->Order - 1)
    return;
  Proposed[DECREASE_ORDER][NewGroup.ID] += 1;
};

void markov::ChangeGroup() {
  group &NewGroup = Groups[Random.irn(0, Groups.size() - 1)];
  if (NewGroup.ID == Var.CurrGroup->ID)
    return;
  if (NewGroup.Order != Var.CurrGroup->Order)
    return;
  Proposed[CHANGE_GROUP][NewGroup.ID] += 1;
};

void markov::ChangeTau(){};

void markov::ChangeMomentum(){};

double markov::RandomPickTau(const double &OldTau, double &Prop) {
  double x = Random.urn();
  double NewTau;
  if (x < 1.0 / 3) {
    double DeltaT = Para.Beta / 3.0;
    NewTau = OldTau + DeltaT * (Random.urn() - 0.5);
  } else if (x < 2.0 / 3) {
    NewTau = -OldTau;
  } else {
    NewTau = Random.urn() * Para.Beta;
  }
  Prop = 1.0;
  if (NewTau < 0.0)
    NewTau += Para.Beta;
  if (NewTau > Para.Beta)
    NewTau -= Para.Beta;
  return NewTau;
};

void markov::RandomPickK(const momentum &OldMom, momentum &NewMom,
                         double &Prop) {
  double x = Random.urn();
  if (x < 1.0 / 3) {
    NewMom = OldMom;
    int dir = Random.irn(0, D - 1);
    double STEP = Para.Beta > 1 ? Para.Kf / Para.Beta * 3.0 : Para.Kf;
    NewMom[dir] += STEP * (Random.urn() - 0.5);
  } else if (x < 2.0 / 3) {
    double k = norm2(OldMom);
    if (k < EPS) {
      Prop = 0.0;
      return;
    }
    const double Lambda = 1.5;
    double knew = k / Lambda + Random.urn() * (Lambda - 1.0 / Lambda) * k;
    double Ratio = knew / k;
    for (int i = 0; i < D; i++)
      NewMom[i] = OldMom[i] * Ratio;
    if (D == 2)
      Prop = 1.0;
    else if (D == 3)
      Prop = Ratio;
  } else {
    for (int i = 0; i < D; i++)
      NewMom[i] = -OldMom[i];
    Prop = 1.0;
  }
};

int markov::RandomPickExtK(const int &OldExtMomBin, double &Prop) {
  Prop = 1.0;
  return Random.irn(0, Para.MaxExtMom - 1);
};
