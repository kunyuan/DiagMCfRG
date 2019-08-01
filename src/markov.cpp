//
//  markov.cpp
//
//  Created by Kun Chen on 1/21/19.
//  Copyright (c) 2019 Kun Chen. All rights reserved.
//
#include "markov.h"
#include "utility/fmt/format.h"
#include "utility/fmt/printf.h"
#include <iostream>

extern parameter Para;
extern RandomFactory Random;

using namespace mc;
using namespace diag;
using namespace std;

#define NAME(x) #x
// #define COPYFROMTO(x, y)                                                       \
//   for (int i = 0; i < D; i++)                                                  \
//     y[i] = x[i];

// x=2*i, then PAIR(x)=2*i+1
// x=2*i+1, then PAIR(x)=2*i
#define PAIR(x) (int(x / 2) * 4 + 1 - x)

markov::markov() : Var(Weight.Var), Groups(Weight.Groups) {
  ///==== initialize Weight ============================//
  Weight.ReadDiagrams();

  //===== initialize updates related variable ==========//

  UpdatesName[INCREASE_ORDER] = NAME(INCREASE_ORDER);
  UpdatesName[DECREASE_ORDER] = NAME(DECREASE_ORDER);
  UpdatesName[CHANGE_GROUP] = NAME(CHANGE_GROUP);
  UpdatesName[CHANGE_MOM] = NAME(CHANGE_MOM);
  UpdatesName[CHANGE_TAU] = NAME(CHANGE_TAU);
  UpdatesName[CHANGE_SCALE] = NAME(CHANGE_SCALE);

  // for(int i=0;i<MCUpdates;i++)
  // UpdatesName[(Updates)i]=NAME((Updates))

  InitialArray(&Accepted[0][0], 1.0e-10, MCUpdates * MaxGroupNum);
  InitialArray(&Proposed[0][0], 1.0e-10, MCUpdates * MaxGroupNum);

  ///==== initialize observable =======================//
  for (auto &g : Groups) {
    Polar[g.ID].fill(1.0e-10);
    PolarStatic[g.ID] = 1.0e-10;
  }
  ///=== Do all kinds of test  =======================//
  Weight.StaticTest();
  Weight.DynamicTest();

  ///==== Set Reweighting factor =====================//
  AdjustGroupReWeight();
};

int markov::DynamicTest() { return Weight.DynamicTest(); }

void markov::PrintDeBugMCInfo() {
  string msg;
  msg = string(80, '=') + "\n";
  msg += "\nMC Counter: " + to_string(Para.Counter) + ":\n";
  msg += "Current Group Info:\n " + ToString(*Var.CurrGroup);

  // msg += string(80, '=') + "\n";
  // msg += "GWeight: \n";
  // for (int i = 0; i < Var.CurrGroup->GNum; i++)
  //   msg += ToString(Var.CurrGroup->Diag[0].G[i]->Weight) + "; ";
  // msg += "\n";

  // msg += "VerWeight: \n";
  // for (int i = 0; i < Var.CurrGroup->Ver4Num; i++) {
  //   msg += ToString(Var.CurrGroup->Diag[0].Ver4[i]->Weight[0]) + ", ";
  //   msg += ToString(Var.CurrGroup->Diag[0].Ver4[i]->Weight[1]) + "; ";
  // }
  // msg += "\n";

  msg += string(80, '=') + "\n";
  msg += "LoopMom: \n";
  for (int d = 0; d < D; d++) {
    for (int i = 0; i < Var.CurrGroup->Order + 3; i++)
      msg += ToString(Var.LoopMom[i][d]) + ", ";
    msg += "\n";
  }
  msg += "\n";

  msg += string(80, '=') + "\n";
  msg += "Tau: \n";
  for (int i = 0; i < Var.CurrGroup->Order + 1; i++)
    msg += ToString(Var.Tau[i]) + ", ";

  msg += "\n";

  LOG_INFO(msg);
}

void markov::AdjustGroupReWeight() {
  for (int i = 0; i < Weight.Groups.size(); i++)
    Weight.Groups[i].ReWeight = Para.ReWeight[i];
};

void markov::Measure() {
  double MCWeight = fabs(Var.CurrGroup->Weight) * Var.CurrGroup->ReWeight;
  double WeightFactor = Var.CurrGroup->Weight / MCWeight;

  Polar[Var.CurrGroup->ID][Var.CurrExtMomBin] += WeightFactor;
  PolarStatic[Var.CurrGroup->ID] += WeightFactor;

  Weight.Measure(1.0 / MCWeight);
  // Weight.Measure(WeightFactor);
};

void markov::UpdateWeight(double Ratio) { Weight.Update(Ratio); }

void markov::SaveToFile() {
  // return;

  for (auto &group : Groups) {
    ofstream PolarFile;
    string FileName = fmt::format("group{0}_pid{1}.dat", group.Name, Para.PID);
    PolarFile.open(FileName, ios::out | ios::trunc);
    if (PolarFile.is_open()) {
      PolarFile << fmt::sprintf(
          "#PID:%d, Type:%d, rs:%.3f, Beta: %.3f, Group: %s, Step: %d\n",
          Para.PID, Para.ObsType, Para.Rs, Para.Beta, group.Name, Para.Counter);

      for (int j = 0; j < Polar[group.ID].size(); j++)
        PolarFile << fmt::sprintf("%13.6f\t%13.6f\n", Para.ExtMomTable[j][0],
                                  Polar[group.ID][j]);
      PolarFile.close();
    } else {
      LOG_WARNING("Polarization for PID " << Para.PID << " fails to save!");
    }
  }

  ofstream StaticPolarFile;
  string FileName = fmt::sprintf("output%d.dat", Para.PID);
  StaticPolarFile.open(FileName, ios::out | ios::trunc);
  if (StaticPolarFile.is_open()) {
    for (auto &group : Groups) {
      StaticPolarFile << fmt::sprintf(
          "PID:%-4d  Type:%-2d  Group:%-4s  rs:%-.3f  "
          "Beta:%-.3f  Lambda:%-.3f  Polar: % 13.6f\n",
          Para.PID, Para.ObsType, group.Name, Para.Rs, Para.Beta, Para.Mass2,
          PolarStatic[group.ID]);
    }
    StaticPolarFile.close();
  } else {
    LOG_WARNING("Static Polarization for PID " << Para.PID
                                               << " fails to save!");
  }

  Weight.Save();
};

void markov::ClearStatis() { Weight.ClearStatis(); }

void markov::ChangeGroup() {
  group &NewGroup = Groups[Random.irn(0, Groups.size() - 1)];
  if (NewGroup.ID == Var.CurrGroup->ID)
    return;

  Updates Name;
  double Prop = 1.0;

  if (NewGroup.Order == Var.CurrGroup->Order) {
    // change group with the same order
    Name = CHANGE_GROUP;
    Prop = 1.0;

  } else if (NewGroup.Order == Var.CurrGroup->Order + 1) {
    // change to a new group with one higher order
    Name = INCREASE_ORDER;
    static momentum NewMom;
    double NewTau1, NewTau2;
    // Generate New Tau
    Prop = GetNewTau(NewTau1, NewTau2);
    int NewTauIndex = Var.CurrGroup->TauNum;
    Var.Tau[NewTauIndex] = NewTau1;
    Var.Tau[NewTauIndex + 1] = NewTau2;
    // Generate New Mom
    Prop *= GetNewK(NewMom);
    Var.LoopMom[Var.CurrGroup->LoopNum] = NewMom;
  } else if (NewGroup.Order == Var.CurrGroup->Order - 1) {
    // change to a new group with one lower order
    Name = DECREASE_ORDER;
    // Remove OldTau
    int TauToRemove = Var.CurrGroup->TauNum - 2;
    Prop = RemoveOldTau(Var.Tau[TauToRemove], Var.Tau[TauToRemove + 1]);
    // Remove OldMom
    int LoopToRemove = Var.CurrGroup->LoopNum - 1;
    Prop *= RemoveOldK(Var.LoopMom[LoopToRemove]);

  } else {
    return;
  }

  Proposed[Name][Var.CurrGroup->ID] += 1;

  // Weight.ChangeGroup(NewGroup);
  double NewWeight = Weight.GetNewWeight(NewGroup) * NewGroup.ReWeight;
  double R = Prop * fabs(NewWeight) / fabs(Var.CurrGroup->Weight) /
             Var.CurrGroup->ReWeight;

  // if (Name == INCREASE_ORDER) {
  //   cout << "time:" << Var.Tau[2] << ", " << Var.Tau[3] << endl;
  //   cout << NewWeight << endl;
  // }
  // if (NewGroup.Order == 2)
  //   cout << NewWeight << endl;

  if (Random.urn() < R) {
    Accepted[Name][Var.CurrGroup->ID]++;
    Weight.AcceptChange(NewGroup);
  } else {
    Weight.RejectChange(NewGroup);
  }
  return;
};

void markov::ChangeTau() {
  int TauIndex = Random.irn(0, Var.CurrGroup->TauNum);

  // if TauIndex is a locked tau, skip
  if (Var.CurrGroup->IsLockedTau[TauIndex]) {
    return;
  }

  Proposed[CHANGE_TAU][Var.CurrGroup->ID]++;

  double CurrTau = Var.Tau[TauIndex];
  double NewTau;
  double Prop = ShiftTau(CurrTau, NewTau);

  Var.Tau[TauIndex] = NewTau;

  // Weight.ChangeTau(*Var.CurrGroup, TauIndex);
  double NewWeight = Weight.GetNewWeight(*Var.CurrGroup);
  double R = Prop * fabs(NewWeight) / fabs(Var.CurrGroup->Weight);
  if (Random.urn() < R) {
    Accepted[CHANGE_TAU][Var.CurrGroup->ID]++;
    Weight.AcceptChange(*Var.CurrGroup);
  } else {
    // retore the old Tau if the update is rejected
    // if TauIndex is external, then its partner can be different
    Var.Tau[TauIndex] = CurrTau;
    Weight.RejectChange(*Var.CurrGroup);
  }
};

void markov::ChangeMomentum() {
  int LoopIndex = Random.irn(0, Var.CurrGroup->LoopNum - 1);
  // int LoopIndex = int(Random.urn() * (Var.CurrGroup->Order + 3));

  // if loopIndex is a locked loop, skip
  if (Var.CurrGroup->IsLockedLoop[LoopIndex])
    return;

  Proposed[CHANGE_MOM][Var.CurrGroup->ID]++;

  double Prop;
  int NewExtMomBin;
  static momentum CurrMom;

  CurrMom = Var.LoopMom[LoopIndex];

  if (Var.CurrGroup->IsExtTransferLoop[LoopIndex]) {
    Prop = ShiftExtTransferK(Var.CurrExtMomBin, NewExtMomBin);
    Var.LoopMom[LoopIndex] = Para.ExtMomTable[NewExtMomBin];
    if (Var.LoopMom[LoopIndex].norm() > Para.MaxExtMom) {
      Var.LoopMom[LoopIndex] = CurrMom;
      return;
    }
  } else if (Var.CurrGroup->IsExtLegLoop[LoopIndex]) {
    Prop = ShiftExtLegK(CurrMom, Var.LoopMom[LoopIndex]);
  } else {
    Prop = ShiftK(CurrMom, Var.LoopMom[LoopIndex]);
  }

  // Weight.ChangeMom(*Var.CurrGroup, LoopIndex);
  double NewWeight = Weight.GetNewWeight(*Var.CurrGroup);
  double R = Prop * fabs(NewWeight) / fabs(Var.CurrGroup->Weight);
  if (Random.urn() < R) {
    Accepted[CHANGE_MOM][Var.CurrGroup->ID]++;
    Weight.AcceptChange(*Var.CurrGroup);
    if (Var.CurrGroup->IsExtTransferLoop[LoopIndex])
      Var.CurrExtMomBin = NewExtMomBin;
  } else {
    Var.LoopMom[LoopIndex] = CurrMom;
    Weight.RejectChange(*Var.CurrGroup);
  }
};

void markov::ChangeScale() {
  if (Para.Type != RG)
    return;
  double Prop = 1.0;
  double OldScale = Var.CurrScale;
  if (Random.urn() < 0.5)
    Var.CurrScale += Random.urn() * 0.5;
  else
    Var.CurrScale -= Random.urn() * 0.5;

  if (Var.CurrScale < 0.0) {
    Var.CurrScale = OldScale;
    return;
  }

  Proposed[CHANGE_SCALE][Var.CurrGroup->ID] += 1;

  // force to change the group weight
  // Weight.ChangeGroup(*(Var.CurrGroup), true);
  double NewWeight = Weight.GetNewWeight(*Var.CurrGroup);

  double R = Prop * fabs(NewWeight) / fabs(Var.CurrGroup->Weight);
  if (Random.urn() < R) {
    Accepted[CHANGE_SCALE][Var.CurrGroup->ID]++;
    Weight.AcceptChange(*Var.CurrGroup);
  } else {
    Var.CurrScale = OldScale;
    Weight.RejectChange(*Var.CurrGroup);
  }
  return;
}

double markov::GetNewTau(double &NewTau1, double &NewTau2) {
  double Step = 1.0;
  NewTau1 = Random.urn() * Para.Beta;
  NewTau2 = Random.urn() * Para.Beta;
  // NewTau2 = (Random.urn() - 0.5) * Step + NewTau1;
  // if (NewTau2 < 0.0)
  //   NewTau2 += Para.Beta;
  // if (NewTau2 > Para.Beta)
  //   NewTau2 -= Para.Beta;
  // return Para.Beta * Step;
  return Para.Beta * Para.Beta;
}

double markov::RemoveOldTau(double &OldTau1, double &OldTau2) {
  // double Step = 1.0;
  // if (abs(OldTau2 - OldTau1) > Step / 2.0)
  //   return 0.0;
  // else
  //   return 1.0 / Para.Beta / Step;
  return 1.0 / Para.Beta / Para.Beta;
}

double markov::GetNewK(momentum &NewMom) {
  //====== The hard Way ======================//
  // double dK = Para.Kf / sqrt(Para.Beta) / 4.0;
  // if (dK > Para.Kf / 2)
  // dK = Para.Kf / 2; // to avoid dK>Kf situation
  // double dK = 1.0 * Var.CurrScale;
  // double x = Random.urn();
  // double KAmp = dK / 2.0 + dK * x;
  // if (x < 0)
  //   KAmp = Para.Kf + KAmp;
  // else
  //   KAmp = Para.Kf - KAmp;
  double dK = Para.Kf / 2.0;
  double KAmp = Para.Kf + (Random.urn() - 0.5) * 2.0 * dK;
  if (KAmp <= 0.0) {
    return 0.0;
  }
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
  // double dK = 1.0 * Var.CurrScale;
  // double dK = Para.Kf / sqrt(Para.Beta) / 4.0;
  // if (dK > Para.Kf / 2)
  //   dK = Para.Kf / 2; // to avoid dK>Kf situation
  double dK = Para.Kf / 2.0;
  double KAmp = OldMom.norm();
  if (KAmp < Para.Kf - dK || KAmp > Para.Kf + dK)
    // if (KAmp < Para.Kf - 1.5 * dK ||
    //     (KAmp > Para.Kf - 0.5 * dK && KAmp < Para.Kf + 0.5 * dK) ||
    //     KAmp > Para.Kf + 1.5 * dK)
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
    if (k < 1.0e-9) {
      Prop = 0.0;
      NewMom = OldMom;
    } else {
      const double Lambda = 1.5;
      double knew = k / Lambda + Random.urn() * (Lambda - 1.0 / Lambda) * k;
      double Ratio = knew / k;
      for (int i = 0; i < D; i++)
        NewMom[i] = OldMom[i] * Ratio;
      if (D == 2)
        Prop = 1.0;
      else if (D == 3)
        Prop = Ratio;
    }

    // if (isnan(Var.LoopMom[0][0]) || isnan(Var.LoopMom[0][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[1][0]) || isnan(Var.LoopMom[1][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[2][0]) || isnan(Var.LoopMom[2][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[3][0]) || isnan(Var.LoopMom[3][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[4][0]) || isnan(Var.LoopMom[4][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[5][0]) || isnan(Var.LoopMom[5][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[6][0]) || isnan(Var.LoopMom[6][0])) {
    //   cout << "Na" << endl;
    // }
  } else {
    for (int i = 0; i < D; i++)
      NewMom[i] = -OldMom[i];
    Prop = 1.0;
  }

  return Prop;
};

double markov::ShiftExtTransferK(const int &OldExtMomBin, int &NewExtMomBin) {
  NewExtMomBin = Random.irn(0, ExtMomBinSize - 1);
  return 1.0;
};

double markov::ShiftExtLegK(const momentum &OldExtMom, momentum &NewExtMom) {
  if (D == 2) {
    double Theta = Random.urn() * 2.0 * PI;
    NewExtMom[0] = Para.Kf * cos(Theta);
    NewExtMom[1] = Para.Kf * sin(Theta);
    return 1.0;
  } else {
    ABORT("ShiftExtKOnKf for D=3 has not yet been implemented!");
  }
};

double markov::ShiftTau(const double &OldTau, double &NewTau) {
  double x = Random.urn();
  if (x < 1.0 / 3) {
    double DeltaT = Para.Beta / 10.0;
    NewTau = OldTau + DeltaT * (Random.urn() - 0.5);
  } else if (x < 2.0 / 3) {
    NewTau = -OldTau;
  } else {
    NewTau = Random.urn() * Para.Beta;
  }

  // NewTau = Random.urn() * Para.Beta;
  if (NewTau < 0.0)
    NewTau += Para.Beta;
  if (NewTau > Para.Beta)
    NewTau -= Para.Beta;
  return 1.0;
}

std::string markov::_DetailBalanceStr(Updates op) {
  string Output = string(80, '-') + "\n";
  Output += UpdatesName[op] + ":\n";
  double TotalProposed = 0.0, TotalAccepted = 0.0;
  for (int i = 0; i <= Groups.size(); i++) {
    if (!Equal(Proposed[op][i], 0.0)) {
      TotalAccepted += Accepted[op][i];
      TotalProposed += Proposed[op][i];
      Output += fmt::sprintf("\t%8s%2i:%15g%15g%15g\n", "Group", Groups[i].ID,
                             Proposed[op][i], Accepted[op][i],
                             Accepted[op][i] / Proposed[op][i]);
      // fmt::format("\t%8s%4s:%15g%15g%15g\n", "Group", Groups[i].Name,
      //             Proposed[op][i], Accepted[op][i],
      //             Accepted[op][i] / Proposed[op][i]);
    }
  }
  if (!Equal(TotalProposed, 0.0)) {
    Output += fmt::sprintf("\t%10s:%15g%15g%15g\n", "Summation", TotalProposed,
                           TotalAccepted, TotalAccepted / TotalProposed);
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
