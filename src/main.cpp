//
//  main.cpp
//
//  Created by Kun Chen on 1/21/19.
//  Copyright (c) 2019 Kun Chen. All rights reserved.
//

/********************** include files *****************************************/
#include "global.h"
#include "markov.h"
#include "utility/abort.h"
#include "utility/logger.h"
#include "utility/timer.h"
#include "weight.h"
#include <iostream>
#include <math.h>
#include <unistd.h>

using namespace std;
using namespace mc;
void InitPara();
void MonteCarlo();

parameter Para; // parameters as a global variable
RandomFactory Random;

int main(int argc, const char *argv[]) {
  cout << "Beta, Rs, Mass2, MaxExtMom(*kF), TotalStep(*1e6), Seed, "
          "PID\n";
  cin >> Para.Beta >> Para.Rs >> Para.Mass2 >> Para.MaxExtMom >>
      Para.TotalStep >> Para.Seed >> Para.PID;
  InitPara(); // initialize global parameters
  MonteCarlo();
  return 0;
}

void InitPara() {
  //// initialize the global log configuration   /////////////
  string LogFile = "_" + to_string(Para.PID) + ".log";
  LOGGER_CONF(LogFile, "MC", Logger::file_on | Logger::screen_on, INFO, INFO);

  // Para.Type = POLAR;
  Para.Type = RG;
  Para.ObsType = FREQ;

  Para.UseVer4 = false;
  // Para.UseVer4 = true;

  // diagram file path: groups/DiagPolar1.dat
  // Para.DiagFileFormat = "groups/DiagPolar{}.txt";
  Para.DiagFileFormat = "groups/DiagLoop{}.txt";
  // Para.GroupName = {"0", "1", "2", "3"};
  // Para.GroupName = {"1", "2"};
  Para.GroupName = {
      "0",
      "1",
      "2",
  };
  Para.ReWeight = {1, 5.0, 3.0, 1.0, 1.0};
  // Para.SelfEnergyType = FOCK;
  Para.SelfEnergyType = BARE;

  Para.Vertex4Type = MOM_ANGLE;

  //// initialize the global parameter //////////////////////
  double Kf;
  if (D == 3) {
    Kf = pow(9.0 * PI / 4.0, 1.0 / 3.0) / Para.Rs; // 3D
  } else if (D == 2) {
    Kf = sqrt(2.0) / Para.Rs; // 2D
  } else {
    ABORT("Dimension " << D << " has not yet been implemented!");
  }
  Para.Kf = Kf;
  Para.Ef = Kf * Kf;
  Para.Mu = Para.Ef;
  Para.MaxExtMom *= Kf;

  // scale all energy with E_F
  Para.Beta /= Para.Ef;
  Para.UVScale = 2.0 * Para.Kf;
  Para.UVCoupling = 1.0 * Para.Ef;

  double dScale = Para.UVScale / ScaleBinSize;
  for (int i = 0; i < (ScaleBinSize + 1); i++) {
    Para.ScaleTable[i] = i * dScale;
    Para.dScaleTable[i] = dScale;
  }
  // double R = 1.5;
  // Para.ScaleTable[0] = 0.0;
  // Para.ScaleTable[1] = 0.4;
  // Para.dScaleTable[0] = Para.ScaleTable[1] - Para.ScaleTable[0];
  // for (int i = 2; i < ScaleBinSize + 1; i++) {
  //   Para.ScaleTable[i] = Para.ScaleTable[i - 1] * R;
  //   Para.dScaleTable[i - 1] = Para.ScaleTable[i] - Para.ScaleTable[i - 1];
  // }

  // double deltaScale = Para.Kf / 16.0;
  // double Start = 0;
  // double End = ScaleBinSize / 2;
  // double dScale = deltaScale / (End - Start);
  // for (int i = Start; i < End; i++) {
  //   Para.ScaleTable[i] = i * dScale;
  // }

  // deltaScale = Para.Kf * 7.0 / 16.0;
  // Start = End;
  // End = ScaleBinSize * 3 / 4;
  // dScale = deltaScale / (End - Start);
  // for (int i = Start; i < End; i++) {
  //   Para.ScaleTable[i] = Para.ScaleTable[i - 1] + dScale;
  // }

  // deltaScale = Para.Kf / 2.0;
  // Start = End;
  // End = ScaleBinSize * 7 / 8;
  // dScale = deltaScale / (End - Start);
  // for (int i = Start; i < End; i++) {
  //   Para.ScaleTable[i] = Para.ScaleTable[i - 1] + dScale;
  // }

  // deltaScale = Para.UVScale - Para.Kf;
  // Start = End;
  // End = ScaleBinSize + 1;
  // dScale = deltaScale / (End - Start);
  // for (int i = Start; i < End; i++) {
  //   Para.ScaleTable[i] = Para.ScaleTable[i - 1] + dScale;
  // }

  // Para.ScaleTable = {0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 8.0};
  // Para.ScaleTable = {0,   0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.2,
  //                    0.3, 0.5,   0.8,  1.0,  1.2,  1.5,  2.0,  8.0};

  Para.ScaleTable[0] = 1.0e-6;
  for (int i = 0; i < ScaleBinSize + 1; i++) {
    Para.dScaleTable[i] = Para.ScaleTable[i + 1] - Para.ScaleTable[i];
  }

  for (int i = 0; i < AngBinSize; i++) {
    Para.AngleTable[i] = diag::Index2Angle(i, AngBinSize);
    Para.dAngleTable[i] = 2.0 * PI / AngBinSize;
  }

  // initialize external momentum
  for (int i = 0; i < ExtMomBinSize; i++) {
    // the external momentum only has x component
    Para.ExtMomTable[i][0] = (i + 0.5) * Para.MaxExtMom / ExtMomBinSize;
    for (int j = 1; j < D; j++)
      Para.ExtMomTable[i][j] = 0.0;
  }

  LOG_INFO("Inverse Temperature: " << Para.Beta << "\n"
                                   << "UV Energy Scale: " << Para.UVScale
                                   << "\n"
                                   << "UV Coupling: " << Para.UVCoupling << "\n"
                                   << "r_s: " << Para.Rs << "\n"
                                   << "Fermi Mom: " << Para.Kf << "\n"
                                   << "Fermi Energy: " << Para.Ef << "\n");

  Para.PrinterTimer = 5;
  Para.SaveFileTimer = 60;
  Para.ReweightTimer = 30;
}

void MonteCarlo() {
  // LOG_INFO("Initializing Markov!");
  markov Markov;
  InterruptHandler Interrupt;

  Random.Reset(Para.Seed);
  Para.Counter = 0;

  timer ReweightTimer, PrinterTimer, SaveFileTimer, MessageTimer;
  PrinterTimer.start();
  SaveFileTimer.start();
  MessageTimer.start();
  ReweightTimer.start();

  LOG_INFO("Start simulation ...")
  long int WaitStep = 1000000;
  int Flag = 0;

  for (int Block = 0; Block < Para.TotalStep; Block++) {
    for (int i = 0; i < 1000000; i++) {
      Para.Counter++;
      // if (Para.Counter == 140737351830544) {
      //   cout << "Before: " << Para.Counter << endl;
      //   Markov.PrintDeBugMCInfo();
      // }

      double x = Random.urn();
      if (x < 1.0 / 4.0) {
        Markov.ChangeGroup();
        // ;
      } else if (x < 2.0 / 4.0) {
        Markov.ChangeMomentum();
        // ;
      } else if (x < 3.0 / 4.0) {
        Markov.ChangeTau();
        // ;
        // } else if (x < 4.0 / 4.0) {
        //   Markov.ChangeScale();
        // ;
      }

      // if (Markov.Weight.Var.CurrScale > Para.UVScale) {
      //   cout << "Scale wrong!" << Markov.Weight.Var.CurrScale << endl;
      //   ABORT("wrong");
      // }

      // if (Para.Counter == 140737351830544) {
      //   cout << "After: " << Para.Counter << endl;
      //   Markov.PrintDeBugMCInfo();
      // }

      if (i % 2 == 0)
        Markov.Measure();
      // Markov.DynamicTest();

      if (i % 1000 == 0) {
        // cout << Markov.Weight.Var.Tau[0] << " vs " <<
        // Markov.Weight.Var.Tau[1]
        //      << endl;
        // cout << Markov.Weight.Var.Tau[2] << " , " << Markov.Weight.Var.Tau[3]
        //      << endl;

        // Markov.PrintDeBugMCInfo();
        if (PrinterTimer.check(Para.PrinterTimer)) {
          Markov.DynamicTest();
          Markov.PrintDeBugMCInfo();
          Markov.PrintMCInfo();
          LOG_INFO(ProgressBar((double)Block / Para.TotalStep));
        }

        if (SaveFileTimer.check(Para.SaveFileTimer)) {
          Interrupt.Delay(); // the process can not be killed in saving
          Markov.SaveToFile();
          Interrupt.Resume(); // after this point, the process can be killed
        }

        if (ReweightTimer.check(Para.ReweightTimer)) {
          Markov.AdjustGroupReWeight();
          Para.ReweightTimer *= 1.5;
        }
      }
    }
    if (Block % 10 == 0) {
      // if (Flag == 0)
      Markov.UpdateWeight(1.0);
      LOG_INFO("Update weight, " << Block);
      // Flag = 1;
      // Markov.ClearStatis();
    }
    // if (i % (WaitStep * 10)) {
    //   // LOG_INFO("Current IR Scale: " << Markov.Var.CurrIRScaleBin);
    // }
  }

  LOG_INFO("Simulation is ended!");
  Markov.PrintMCInfo();
  Interrupt.Delay(); // the process can not be killed in saving
  Markov.SaveToFile();
  Interrupt.Resume(); // after this point, the process can be killed
  LOG_INFO("Quit Markov.");
}
