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
#define COPYFROMTO(x, y)                                                       \
  for (int i; i < D; i++)                                                      \
    y[i] = x[i];

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
  Proposed[INCREASE_ORDER][Var.CurrGroup->ID] += 1;

  // Generate New Tau
  double NewTau;
  double Prop = GetNewTau(NewTau);
  int NewTauIndex = Var.CurrGroup->TauNum;
  Var.Tau[NewTauIndex] = NewTau;
  Var.Tau[NewTauIndex + 1] = NewTau;

  // Generate New Mom
  static momentum NewMom;
  Prop *= GetNewK(NewMom);
  int NewLoopIndex = Var.CurrGroup->LoopNum;
  COPYFROMTO(NewMom, Var.LoopMom[NewLoopIndex]);

  Weight.ChangeGroup(NewGroup);
  double NewWeight = Weight.GetNewWeight(NewGroup);
  double R = Prop * fabs(NewWeight) / fabs(Var.CurrWeight);
  if (Random.urn() < R) {
    Accepted[INCREASE_ORDER][Var.CurrGroup->ID]++;
    Weight.AcceptChange(NewGroup);
  }
};

void markov::DecreaseOrder() {
  group &NewGroup = Groups[Random.irn(0, Groups.size() - 1)];
  if (NewGroup.Order != Var.CurrGroup->Order - 1)
    return;
  Proposed[DECREASE_ORDER][Var.CurrGroup->ID] += 1;

  // Remove OldTau
  int TauToRemove = Var.CurrGroup->TauNum - 1;
  double Prop = RemoveOldTau(Var.Tau[TauToRemove]);
  // Remove OldMom
  int LoopToRemove = Var.CurrGroup->LoopNum - 1;
  Prop *= RemoveOldK(Var.LoopMom[LoopToRemove]);

  Weight.ChangeGroup(NewGroup);
  double NewWeight = Weight.GetNewWeight(NewGroup);
  double R = Prop * fabs(NewWeight) / fabs(Var.CurrWeight);
  if (Random.urn() < R) {
    Accepted[DECREASE_ORDER][Var.CurrGroup->ID]++;
    Weight.AcceptChange(NewGroup);
  }
};

void markov::ChangeGroup() {
  group &NewGroup = Groups[Random.irn(0, Groups.size() - 1)];
  if (NewGroup.ID == Var.CurrGroup->ID)
    return;
  if (NewGroup.Order != Var.CurrGroup->Order)
    return;
  Proposed[INCREASE_ORDER][Var.CurrGroup->ID] += 1;
  return;
};

void markov::ChangeTau() {
  // TauIndex can be 0, 2, 4, ...
  int TauIndex = Random.irn(0, Var.CurrGroup->TauNum / 2 - 1) * 2;
  if (Para.ObsType == 1 && TauIndex == 0)
    return;
  Proposed[CHANGE_TAU][Var.CurrGroup->ID]++;

  // note if TauIndex==0, then Tau[1]!=Tau[0].
  // since we only change Tau[1] in this case, it is better set CurrTau=Tau[Odd]
  double CurrTau = Var.Tau[TauIndex + 1];
  double NewTau;
  double Prop = ShiftTau(CurrTau, NewTau);
  if (TauIndex != 0)
    Var.Tau[TauIndex] = NewTau;
  Var.Tau[TauIndex + 1] = NewTau;

  Weight.ChangeTau(*Var.CurrGroup, TauIndex);
  double NewWeight = Weight.GetNewWeight(*Var.CurrGroup);
  double R = Prop * fabs(NewWeight) / fabs(Var.CurrWeight);
  if (Random.urn() < R) {
    Accepted[CHANGE_TAU][Var.CurrGroup->ID]++;
    Weight.AcceptChange(*Var.CurrGroup);
  } else {
    // retore the old Tau if the update is rejected
    if (TauIndex != 1)
      Var.Tau[TauIndex] = CurrTau;
    Var.Tau[TauIndex + 1] = CurrTau;
  }
};

void markov::ChangeMomentum() {
  static momentum NewMom, CurrMom;
  int LoopIndex = Random.irn(0, Var.CurrGroup->LoopNum - 1);
  Proposed[CHANGE_MOM][Var.CurrGroup->ID]++;

  COPYFROMTO(Var.LoopMom[LoopIndex], CurrMom);

  double Prop;
  if (LoopIndex == 0) {
    int NewExtMomBin;
    Prop = ShiftExtK(Var.CurrExtMomBin, NewExtMomBin);
    COPYFROMTO(Var.ExtMomTable[NewExtMomBin], NewMom);
  } else {
    Prop = ShiftK(CurrMom, NewMom);
  }
  if (LoopIndex == 0 && norm2(NewMom) > Para.MaxExtMom)
    return;

  Weight.ChangeMom(*Var.CurrGroup, LoopIndex);
  double NewWeight = Weight.GetNewWeight(*Var.CurrGroup);
  double R = Prop * fabs(NewWeight) / fabs(Var.CurrWeight);
  if (Random.urn() < R) {
    Accepted[CHANGE_MOM][Var.CurrGroup->ID]++;
    Weight.AcceptChange(*Var.CurrGroup);
  } else {
    COPYFROMTO(CurrMom, Var.LoopMom[LoopIndex]);
  }
};

double markov::GetNewTau(double &NewTau) {
  NewTau = Random.urn() * Para.Beta;
  return Para.Beta;
};
double markov::RemoveOldTau(double &OldTau) { return 1.0 / Para.Beta; }

double markov::GetNewK(momentum &NewMom) {
  //====== The hard Way ======================//
  double dK = Para.Kf / sqrt(Para.Beta) / 4.0;
  if (dK > Para.Kf / 2)
    dK = Para.Kf / 2; // to avoid dK>Kf situation
  double KAmp = Para.Kf + (Random.urn() - 0.5) * 2.0 * dK;
  // Kf-dK<KAmp<Kf+dK
  double Phi = 2.0 * PI * Random.urn();
  if (D == 3) {
    double Theta = PI * Random.urn();
    if (Theta == 0.0)
      return 0.0;
    double K_XY = KAmp * sin(Theta);
    NewMom[0] = K_XY * cos(Phi);
    NewMom[1] = K_XY * sin(Phi);
    NewMom[D - 1] = K_XY * cos(Theta);
    return 2.0 * dK                    // prop density of KAmp in [Kf-dK, Kf+dK)
           * 2.0 * PI                  // prop density of Phi
           * PI                        // prop density of Theta
           * sin(Theta) * KAmp * KAmp; // Jacobian
  } else if (D == 2) {
    NewMom[0] = KAmp * cos(Phi);
    NewMom[1] = KAmp * sin(Phi);
    return 2.0 * dK   // prop density of KAmp in [Kf-dK, Kf+dK)
           * 2.0 * PI // prop density of Phi
           * KAmp;    // Jacobian
  }

  //===== The simple way  =======================//
  // for (int i = 0; i < D; i++)
  //   NewMom[i] = Para.Kf * (Random.urn() - 0.5) * 2.0;
  // return pow(2.0 * Para.Kf, D);

  //============================================//
};

double markov::RemoveOldK(momentum &OldMom) {
  //====== The hard Way ======================//
  double dK = Para.Kf / sqrt(Para.Beta) / 4.0;
  if (dK > Para.Kf / 2)
    dK = Para.Kf / 2; // to avoid dK>Kf situation
  double KAmp = norm2(OldMom);
  if (KAmp < Para.Kf - dK || KAmp > Para.Kf + dK)
    // Kf-dK<KAmp<Kf+dK
    return 0.0;
  if (D == 3) {
    auto SinTheta = sqrt(OldMom[0] * OldMom[0] + OldMom[1] * OldMom[1]) / KAmp;
    if (SinTheta < EPS)
      return 0.0;
    return 1.0 / (2.0 * dK * 2.0 * PI * PI * SinTheta * KAmp * KAmp);
  } else if (D == 2) {
    return 1.0 / (2.0 * dK * 2.0 * PI * KAmp);
  }

  //===== The simple way  =======================//
  // for (int i = 0; i < D; i++)
  //   if(fabs(OldMom[i]>Para.Kf)
  //        return 0.0;
  // return 1.0/pow(2.0 * Para.Kf, D);
  //============================================//
}

double markov::ShiftK(const momentum &OldMom, momentum &NewMom) {
  double x = Random.urn();
  double Prop;
  if (x < 1.0 / 3) {
    NewMom = OldMom;
    int dir = Random.irn(0, D - 1);
    double STEP = Para.Beta > 1 ? Para.Kf / Para.Beta * 3.0 : Para.Kf;
    NewMom[dir] += STEP * (Random.urn() - 0.5);

  } else if (x < 2.0 / 3) {
    double k = norm2(OldMom);
    if (k < EPS)
      Prop = 0.0;

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
  return Prop;
};

double markov::ShiftExtK(const int &OldExtMomBin, int &NewExtMomBin) {
  NewExtMomBin = Random.irn(0, Para.MaxExtMom - 1);
  return 1.0;
};

double markov::ShiftTau(const double &OldTau, double &NewTau) {
  double x = Random.urn();
  if (x < 1.0 / 3) {
    double DeltaT = Para.Beta / 3.0;
    NewTau = OldTau + DeltaT * (Random.urn() - 0.5);
  } else if (x < 2.0 / 3) {
    NewTau = -OldTau;
  } else {
    NewTau = Random.urn() * Para.Beta;
  }
  if (NewTau < 0.0)
    NewTau += Para.Beta;
  if (NewTau > Para.Beta)
    NewTau -= Para.Beta;
  return 1.0;
};