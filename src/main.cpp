//
//  main.cpp
//
//  Created by Kun Chen on 1/21/19.
//  Copyright (c) 2019 Kun Chen. All rights reserved.
//

/********************** include files *****************************************/
#include <iostream>
#include <unistd.h>
#include <math.h>
#include "utility/timer.h"
#include "utility/logger.h"
#include "utility/abort.h"
#include "global.h"
#include "weight.h"
#include "markov.h"

using namespace std;
int RunTest();
void Initilization();
void MonteCarlo();

parameter Para; // parameters as a global variable
RandomFactory Random;

int main(int argc, const char *argv[])
{
    cout << "Beta, Rs, Mass2, MaxExtMom(*kF), TotalStep(*1e6), Observable, Seed, PID\n";
    cin >> Para.Beta >> Para.Rs >> Para.Mass2 >> Para.MaxExtMom >> Para.TotalStep >> Para.ObsType >> Para.Seed >> Para.PID;
    Initilization();
    diag::weight Weight; //basis to calculate diagram weight
    Weight.ReadDiagrams("DiagPolar");
    Weight.Initialization(); //initialize MC variables
    RunTest();
    MonteCarlo();
    return 0;
}

void Initilization()
{
    //// initialize the global log configuration   /////////////
    string LogFile = "_" + to_string(Para.PID) + ".log";
    LOGGER_CONF(LogFile, "MC", Logger::file_on | Logger::screen_on, INFO, INFO);

    Random.Reset(Para.Seed);

    //// initialize the global parameter //////////////////////
    double Kf;
    if (D == 3)
    {
        Kf = pow(9.0 * PI / 4.0, 1.0 / 3.0) / Para.Rs; //3D
    }
    else if (D == 2)
    {
        Kf = sqrt(2.0) / Para.Rs; //2D
    }
    else
    {
        ABORT("Dimension " << D << " has not yet been implemented!");
    }
    Para.Kf = Kf;
    Para.Ef = Kf * Kf;
    Para.Mu = Para.Ef;
    Para.MaxExtMom *= Kf;

    Para.UVScale = 8.0 * Para.Ef;
    Para.UVCoupling = 1.0 * Para.Ef;
    LOG_INFO("Inverse Temperature: " << Para.Beta << "\n"
                                     << "UV Energy Scale: " << Para.UVScale << "\n"
                                     << "UV Coupling: " << Para.UVCoupling << "\n"
                                     << "r_s: " << Para.Rs << "\n"
                                     << "Fermi Mom: " << Para.Kf << "\n"
                                     << "Fermi Energy: " << Para.Ef << "\n");
}

#define TEST(func)                  \
    {                               \
        if (EXIT_SUCCESS != func()) \
            exit(0);                \
    }

int RunTest()
{
    //TestTimer();  //Test the timer
    //TestRNG();
    //TestArray();
    //    TEST(mc::TestMarkov);
    //    TEST(mc::TestDiagCounter);

    //    TEST(TestDictionary);

    return 0;
}

void MonteCarlo()
{
    InterruptHandler Interrupt;
    // auto& Markov = Env.Markov;
    // auto& MarkovMonitor = Env.MarkovMonitor;
    // auto& Para = Env.Para;

    LOG_INFO("Markov is started!");
    timer ReweightTimer, PrinterTimer, DiskWriterTimer, MessageTimer;
    PrinterTimer.start();
    DiskWriterTimer.start();
    MessageTimer.start();
    ReweightTimer.start();

    uint Step = 0;
    // while (true) {
    //Don't use Para.Counter as counter
    // Step++;
    // Markov.Hop(Para.Sweep);
    // MarkovMonitor.Measure();
    //        if (!Env.Diag.CheckDiagram()){
    //            ABORT("Diagram Check didn't pass!!!");
    //        }
    // if (!Markov.Diag->Worm.Exist) {
    //     if (Markov.Diag->MeasureGLine)
    //         sigma[Markov.Diag->Order]++;
    //     else
    //         polar[Markov.Diag->Order]++;
    // }

    // if (Step % 100 == 0) {
    //     // MarkovMonitor.AddStatistics();
    //     if (PrinterTimer.check(Para.PrinterTimer)) {
    //         if (!Env.Diag.CheckDiagram()){
    //             ABORT("Diagram Check didn't pass!!!");
    //         }
    //         Markov.PrintDetailBalanceInfo();
    //         MarkovMonitor.PrintOrderReWeight();
    //     }

    //     if (DiskWriterTimer.check(Para.DiskWriterTimer)) {
    //         Interrupt.Delay();
    //         Env.Save();
    //         Interrupt.Resume();
    //     }

    //     if (MessageTimer.check(Para.MessageTimer))
    //         Env.ListenToMessage();

    //     if (ReweightTimer.check(Para.ReweightTimer))
    //         Env.AdjustOrderReWeight();
    // }
    // }
    LOG_INFO("Markov is ended!");
}