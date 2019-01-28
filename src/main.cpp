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
void MonteCarlo(markov &);

parameter Para; // parameters as a global variable
RandomFactory Random;

int main(int argc, const char *argv[]) {
  cout << "Beta, Rs, Mass2, MaxExtMom(*kF), TotalStep(*1e6), Observable, Seed, "
          "PID\n";
  cin >> Para.Beta >> Para.Rs >> Para.Mass2 >> Para.MaxExtMom >>
      Para.TotalStep >> Para.ObsType >> Para.Seed >> Para.PID;
  InitPara();
  markov Markov;
  Markov.Initialization("DiagPolar"); // initialize MC
  MonteCarlo(Markov);
  return 0;
}

void InitPara() {
  //// initialize the global log configuration   /////////////
  string LogFile = "_" + to_string(Para.PID) + ".log";
  LOGGER_CONF(LogFile, "MC", Logger::file_on | Logger::screen_on, INFO, INFO);

  Random.Reset(Para.Seed);

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
  Para.UVScale = 8.0 * Para.Ef;
  Para.UVCoupling = 1.0 * Para.Ef;

  LOG_INFO("Inverse Temperature: " << Para.Beta << "\n"
                                   << "UV Energy Scale: " << Para.UVScale
                                   << "\n"
                                   << "UV Coupling: " << Para.UVCoupling << "\n"
                                   << "r_s: " << Para.Rs << "\n"
                                   << "Fermi Mom: " << Para.Kf << "\n"
                                   << "Fermi Energy: " << Para.Ef << "\n");

  Para.PrinterTimer = 5;
  Para.SaveFileTimer = 30;
  Para.ReweightTimer = 30;

  Para.GroupID = {
      1,
      2,
      3,
  };
}

void MonteCarlo(markov &Markov) {
  InterruptHandler Interrupt;

  LOG_INFO("Markov is started!");
  timer ReweightTimer, PrinterTimer, SaveFileTimer, MessageTimer;
  PrinterTimer.start();
  SaveFileTimer.start();
  MessageTimer.start();
  ReweightTimer.start();

  const int SWEEP = 1;
  int Step = 0;
  while (Step < Para.TotalStep) {
    Step++;
    for (int i = 0; i < 1000000; i++) {
      Para.Counter++;
      // if (Para.Counter == 9) {
      //   cout << "Before: " << Para.Counter << endl;
      //   PrintDeBugMCInfo();
      // }

      double x = Random.urn();
      if (x < 1.0 / MCUpdates)
        Markov.IncreaseOrder();
      else if (x < 2.0 / MCUpdates)
        Markov.DecreaseOrder();
      else if (x < 3.0 / MCUpdates)
        Markov.ChangeGroup();
      else if (x < 4.0 / MCUpdates)
        Markov.ChangeMomentum();
      // ;
      else if (x < 5.0 / MCUpdates)
        Markov.ChangeTau();
      // ;

      // if (Para.Counter == 8831001) {
      //   cout << "After: " << Para.Counter << endl;
      //   PrintDeBugMCInfo();
      // }

      // double Tau = Var.Tau[1] - Var.Tau[0];
      // momentum G1, G2;
      // for (int i = 0; i < D; i++) {
      //   G1[i] = Var.LoopMom[0][i] + Var.LoopMom[1][i];
      //   G2[i] = Var.LoopMom[1][i];
      // }
      // ASSERT_ALLWAYS(Equal(Green(Tau, G1, UP, 0),
      //                      Var.CurrGroup->Diag[0].G[0]->Weight, 1.0e-8),
      //                Weight._DebugInfo());

      // ASSERT_ALLWAYS(Equal(Green(-Tau, G2, UP, 0),
      //                      Var.CurrGroup->Diag[0].G[1]->Weight, 1.0e-8),
      //                Weight._DebugInfo());
      Markov.Measure();

      if (i % 100 == 0) {
        // Markov.PrintDeBugMCInfo();
        if (PrinterTimer.check(Para.PrinterTimer)) {
          Markov.DynamicTest();
          Markov.PrintDeBugMCInfo();
          Markov.PrintMCInfo();
          LOG_INFO(ProgressBar((double)Step / Para.TotalStep));
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
  }

  Markov.PrintMCInfo();
  Interrupt.Delay(); // the process can not be killed in saving
  Markov.SaveToFile();
  Interrupt.Resume(); // after this point, the process can be killed
  LOG_INFO("Markov is ended!");
}