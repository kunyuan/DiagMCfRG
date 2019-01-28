//
//  markov.cpp
//
//  Created by Kun Chen on 1/21/19.
//  Copyright (c) 2019 Kun Chen. All rights reserved.
//
#include "markov.h"
#include "utility/tinyformat.h"
#include <stdio.h>

extern parameter Para;
extern RandomFactory Random;

using namespace mc;
using namespace diag;
using namespace std;

#define NAME(x) #x
#define COPYFROMTO(x, y)                                                       \
  for (int i = 0; i < D; i++)                                                  \
    y[i] = x[i];

int markov::DynamicTest() { return Weight.DynamicTest(); }

void markov::PrintDeBugMCInfo() {
  string msg;
  msg = string(80, '=') + "\n";
  msg += "\nMC Counter: " + to_string(Para.Counter) + ":\n";
  msg += "Current Group Info:\n " + ToString(*Var.CurrGroup);

  msg += string(80, '=') + "\n";
  msg += "GWeight: \n";
  for (int i = 0; i < Var.CurrGroup->GNum; i++)
    msg += ToString(Var.CurrGroup->Diag[0].G[i]->Weight) + "; ";
  msg += "\n";

  msg += "VerWeight: \n";
  for (int i = 0; i < Var.CurrGroup->Ver4Num; i++) {
    msg += ToString(Var.CurrGroup->Diag[0].Ver[i]->Weight[0]) + ", ";
    msg += ToString(Var.CurrGroup->Diag[0].Ver[i]->Weight[1]) + "; ";
  }
  msg += "\n";

  msg += string(80, '=') + "\n";
  msg += "LoopMom: \n";
  for (int d = 0; d < D; d++) {
    for (int i = 0; i < Var.CurrGroup->LoopNum; i++)
      msg += ToString(Var.LoopMom[i][d]) + ", ";
    msg += "\n";
  }
  msg += "\n";

  msg += string(80, '=') + "\n";
  msg += "Tau: \n";
  for (int i = 0; i < Var.CurrGroup->TauNum; i++)
    msg += ToString(Var.Tau[i]) + ", ";

  msg += "\n";

  LOG_INFO(msg);
}

void markov::Initialization(string FilePrefix) {
  ///==== initialize Weight ============================//
  Weight.ReadDiagrams(FilePrefix);
  Weight.Initialization();

  //===== initialize updates related variable ==========//
  Para.Counter = 0;

  UpdatesName[INCREASE_ORDER] = NAME(INCREASE_ORDER);
  UpdatesName[DECREASE_ORDER] = NAME(DECREASE_ORDER);
  UpdatesName[CHANGE_GROUP] = NAME(CHANGE_GROUP);
  UpdatesName[CHANGE_MOM] = NAME(CHANGE_MOM);
  UpdatesName[CHANGE_TAU] = NAME(CHANGE_TAU);

  // for(int i=0;i<MCUpdates;i++)
  // UpdatesName[(Updates)i]=NAME((Updates))

  InitialArray(&Accepted[0][0], 1.0e-10, MCUpdates * MaxGroupNum);
  InitialArray(&Proposed[0][0], 1.0e-10, MCUpdates * MaxGroupNum);

  ///==== initialize observable =======================//
  for (auto &g : Weight.Groups) {
    Polar[g.ID].fill(1.0e-10);
    PolarStatic[g.ID] = 1.0e-10;
  }
  ///=== Do all kinds of test  =======================//
  Weight.StaticTest();
  Weight.DynamicTest();

  ///==== Set Reweighting factor =====================//
  AdjustGroupReWeight();
}

void markov::AdjustGroupReWeight() {
  double factor[3] = {1.0, 1.0, 10.0};
  for (int i = 0; i < Weight.Groups.size(); i++)
    Weight.Groups[i].ReWeight = factor[i];
};

void markov::Measure() {
  double MCWeight = fabs(Var.CurrGroup->Weight) * Var.CurrGroup->ReWeight;
  double WeightFactor = Var.CurrGroup->Weight / MCWeight;

  Polar[Var.CurrGroup->ID][Var.CurrExtMomBin] += WeightFactor;
  PolarStatic[Var.CurrGroup->ID] += WeightFactor;
};

void markov::SaveToFile() {
  for (auto &id : Para.GroupID) {
    ofstream PolarFile;
    string FileName = tfm::format("group%d_pid%d.dat", id, Para.PID);
    PolarFile.open(FileName, ios::out | ios::trunc);
    if (PolarFile.is_open()) {
      PolarFile << tfm::format(
          "#PID: %d, Type:%d, rs:%.3f, Beta: %.3f, Group: %d, Step: %d\n",
          Para.PID, Para.ObsType, Para.Rs, Para.Beta, id, Para.Counter);

      for (int j = 0; j < Polar[id].size(); j++)
        PolarFile << tfm::format("%13.6f\t%13.6f\n", Var.ExtMomTable[j][0],
                                 Polar[id][j]);
      PolarFile.close();
    } else {
      LOG_WARNING("Polarization for PID " << Para.PID << " fails to save!");
    }
  }

  ofstream StaticPolarFile;
  string FileName = tfm::format("output%d.dat", Para.PID);
  StaticPolarFile.open(FileName, ios::out | ios::trunc);
  if (StaticPolarFile.is_open()) {
    for (auto &id : Para.GroupID) {
      StaticPolarFile << tfm::format(
          "PID:%-4d  Type:%-2d  Group:%-4d  rs:%-.3f  "
          "Beta:%-.3f  Lambda:%-.3f  Polar: % 13.6f\n",
          Para.PID, Para.ObsType, id, Para.Rs, Para.Beta, Para.Mass2,
          PolarStatic[id]);
    }
    StaticPolarFile.close();
  } else {
    LOG_WARNING("Static Polarization for PID " << Para.PID
                                               << " fails to save!");
  }
};

void markov::IncreaseOrder() {
  static momentum NewMom;

  group &NewGroup = Groups[Random.irn(0, Groups.size() - 1)];
  if (NewGroup.Order != Var.CurrGroup->Order + 1)
    return;
  Proposed[INCREASE_ORDER][Var.CurrGroup->ID] += 1;

  // Generate New Tau
  double NewTau;
  double Prop = GetNewTau(NewTau);
  int NewTauIndex = Var.CurrGroup->TauNum;
  // ASSUME: NewTauIndex will never equal to 0 or 1
  Var.Tau[NewTauIndex] = NewTau;
  Var.Tau[NewTauIndex + 1] = NewTau;

  // Generate New Mom
  Prop *= GetNewK(NewMom);
  int NewLoopIndex = Var.CurrGroup->LoopNum;
  COPYFROMTO(NewMom, Var.LoopMom[NewLoopIndex]);

  Weight.ChangeGroup(NewGroup);
  double NewWeight = Weight.GetNewWeight(NewGroup) * NewGroup.ReWeight;
  double R = Prop * fabs(NewWeight) / fabs(Var.CurrGroup->Weight) /
             Var.CurrGroup->ReWeight;

  if (Random.urn() < R) {
    Accepted[INCREASE_ORDER][Var.CurrGroup->ID]++;
    Weight.AcceptChange(NewGroup);
  } else {
    Weight.RejectChange(NewGroup);
  }
};

void markov::DecreaseOrder() {
  group &NewGroup = Groups[Random.irn(0, Groups.size() - 1)];
  if (NewGroup.Order != Var.CurrGroup->Order - 1)
    return;
  Proposed[DECREASE_ORDER][Var.CurrGroup->ID] += 1;

  // Remove OldTau
  int TauToRemove = Var.CurrGroup->TauNum - 2;
  double Prop = RemoveOldTau(Var.Tau[TauToRemove]);
  // Remove OldMom
  int LoopToRemove = Var.CurrGroup->LoopNum - 1;
  Prop *= RemoveOldK(Var.LoopMom[LoopToRemove]);

  Weight.ChangeGroup(NewGroup);
  double NewWeight = Weight.GetNewWeight(NewGroup) * NewGroup.ReWeight;
  double R = Prop * fabs(NewWeight) / fabs(Var.CurrGroup->Weight) /
             Var.CurrGroup->ReWeight;

  if (Random.urn() < R) {
    Accepted[DECREASE_ORDER][Var.CurrGroup->ID]++;
    Weight.AcceptChange(NewGroup);
  } else {
    Weight.RejectChange(NewGroup);
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
  // since we only change Tau[1] in this case, it is better set
  // CurrTau=Tau[Odd]
  double CurrTau = Var.Tau[TauIndex + 1];
  double NewTau;
  double Prop = ShiftTau(CurrTau, NewTau);
  if (TauIndex != 0)
    Var.Tau[TauIndex] = NewTau;
  Var.Tau[TauIndex + 1] = NewTau;

  Weight.ChangeTau(*Var.CurrGroup, TauIndex);
  double NewWeight = Weight.GetNewWeight(*Var.CurrGroup);
  double R = Prop * fabs(NewWeight) / fabs(Var.CurrGroup->Weight);
  if (Random.urn() < R) {
    Accepted[CHANGE_TAU][Var.CurrGroup->ID]++;
    Weight.AcceptChange(*Var.CurrGroup);
  } else {
    // retore the old Tau if the update is rejected
    if (TauIndex != 0)
      Var.Tau[TauIndex] = CurrTau;
    Var.Tau[TauIndex + 1] = CurrTau;
    Weight.RejectChange(*Var.CurrGroup);
  }
};

void markov::ChangeMomentum() {
  int LoopIndex = Random.irn(0, Var.CurrGroup->LoopNum - 1);
  Proposed[CHANGE_MOM][Var.CurrGroup->ID]++;

  double Prop;
  int NewExtMomBin;
  static momentum CurrMom;

  // COPYFROMTO(Var.LoopMom[LoopIndex], CurrMom);
  CurrMom = Var.LoopMom[LoopIndex];

  if (LoopIndex == 0) {
    Prop = ShiftExtK(Var.CurrExtMomBin, NewExtMomBin);
    // COPYFROMTO(Var.ExtMomTable[NewExtMomBin], Var.LoopMom[LoopIndex]);
    Var.LoopMom[LoopIndex] = Var.ExtMomTable[NewExtMomBin];
  } else {
    Prop = ShiftK(CurrMom, Var.LoopMom[LoopIndex]);
  }
  if (LoopIndex == 0 && Var.LoopMom[LoopIndex].norm() > Para.MaxExtMom) {
    Var.LoopMom[LoopIndex] = CurrMom;
    // COPYFROMTO(CurrMom, Var.LoopMom[LoopIndex]);
    return;
  }

  Weight.ChangeMom(*Var.CurrGroup, LoopIndex);
  double NewWeight = Weight.GetNewWeight(*Var.CurrGroup);
  double R = Prop * fabs(NewWeight) / fabs(Var.CurrGroup->Weight);
  if (Random.urn() < R) {
    Accepted[CHANGE_MOM][Var.CurrGroup->ID]++;
    Weight.AcceptChange(*Var.CurrGroup);
    if (LoopIndex == 0)
      Var.CurrExtMomBin = NewExtMomBin;
  } else {
    Var.LoopMom[LoopIndex] = CurrMom;
    // COPYFROMTO(CurrMom, Var.LoopMom[LoopIndex]);
    Weight.RejectChange(*Var.CurrGroup);
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
    NewMom[D - 1] = KAmp * cos(Theta);
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
  double KAmp = OldMom.norm();
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
  //   if (fabs(OldMom[i] > Para.Kf))
  //     return 0.0;
  // return 1.0 / pow(2.0 * Para.Kf, D);
  //============================================//
}

double markov::ShiftK(const momentum &OldMom, momentum &NewMom) {
  double x = Random.urn();
  double Prop;
  if (x < 1.0 / 3) {
    // COPYFROMTO(OldMom, NewMom);
    NewMom = OldMom;
    int dir = Random.irn(0, D - 1);
    double STEP = Para.Beta > 1.0 ? Para.Kf / Para.Beta * 3.0 : Para.Kf;
    NewMom[dir] += STEP * (Random.urn() - 0.5);
    Prop = 1.0;
  } else if (x < 2.0 / 3) {
    double k = OldMom.norm();
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
  NewExtMomBin = Random.irn(0, ExtMomBinSize - 1);
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
}

std::string markov::_DetailBalanceStr(Updates op) {
  string Output = string(80, '-') + "\n";
  Output += UpdatesName[op] + ":\n";
  char temp[80];
  double TotalProposed = 0.0, TotalAccepted = 0.0;
  for (int i = 0; i <= Groups.size(); i++) {
    if (!Equal(Proposed[op][i], 0.0)) {
      TotalAccepted += Accepted[op][i];
      TotalProposed += Proposed[op][i];
      sprintf(temp, "\t%8s%2i:%15g%15g%15g\n", "Group", i, Proposed[op][i],
              Accepted[op][i], Accepted[op][i] / Proposed[op][i]);
      Output += temp;
    }
  }
  if (!Equal(TotalProposed, 0.0)) {
    sprintf(temp, "\t%10s:%15g%15g%15g\n", "Summation", TotalProposed,
            TotalAccepted, TotalAccepted / TotalProposed);
    Output += temp;
  } else
    Output += "\tNo updates are proposed/accepted!\n";
  return Output;
}

void markov::PrintMCInfo() {
  string Output = "";
  Output = string(80, '=') + "\n";
  Output += "MC Counter: " + to_string(Para.Counter) + "\n";
  for (int i = 0; i < MCUpdates; i++)
    Output += _DetailBalanceStr((Updates)i);
  Output += string(80, '=') + "\n";
  LOG_INFO(Output);
}
