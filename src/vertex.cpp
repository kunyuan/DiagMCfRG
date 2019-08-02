#include "vertex.h"
#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/fmt/printf.h"
#include "utility/utility.h"
#include <cmath>
#include <iostream>

using namespace diag;
using namespace std;

extern parameter Para;

// double sum2(const momentum &Mom) {
//   double Sum2 = 0.0;
//   for (int i = 0; i < D; i++)
//     Sum2 += Mom[i] * Mom[i];
//   return Sum2;
// }

double verQ::Interaction(double Tau, const momentum &Mom, int VerType,
                         double Scale) {
  if (VerType >= 0) {
    double interaction = 8.0 * PI / (Mom.squaredNorm() + Para.Mass2);
    if (VerType > 0) {
      // the interaction contains counter-terms
      interaction *=
          pow(Para.Mass2 / (Mom.squaredNorm() + Para.Mass2), VerType);
      interaction *= pow(-1, VerType);
      return interaction;
    }
  } else if (VerType == -1) {
    return 1.0;
  } else if (VerType == -2) {
    return 0.0;
  } else {
    ABORT("VerType can not be " << VerType);
  }
}

// double norm2(const momentum &Mom) { return sqrt(sum2(Mom)); }
verQTheta::verQTheta() {

  if (D == 3)
    return;

  _TestAngle2D();
  _TestAngleIndex();

  Normalization = 1.0e-10;
  // PhyWeight = Para.Kf * (1.0 - exp(-Para.MaxExtMom / Para.Kf)) * 4.0 * PI *
  // PI;
  // PhyWeight = (1.0 - exp(-Para.MaxExtMom / Para.Kf)) /
  //             (1.0 - exp(-Para.MaxExtMom / Para.Kf / ExtMomBinSize)) * 4.0 *
  //             PI * PI;
  // PhyWeight =
  //     1.0 / Para.Beta / Para.Beta * ExtMomBinSize * 2.0 * PI * Para.Kf * 4.0;

  PhyWeight = ExtMomBinSize * 2.0 * PI * Para.Beta;
  // PhyWeight = 2.0 * PI / TauBinSize * 64;
  // PhyWeight = 1.0;

  QIndex = TauBinSize;
  AngleIndex = TauBinSize * ExtMomBinSize;
  OrderIndex = TauBinSize * ExtMomBinSize * AngBinSize;

  EffInteraction = new double[AngBinSize * ExtMomBinSize * TauBinSize];

  DiffInteraction =
      new double[MaxOrder * AngBinSize * ExtMomBinSize * TauBinSize];
  // double DiffInteraction[MaxOrder][ScaleBinSize +
  // 1][AngBinSize][ExtMomBinSize]; double IntInteraction[MaxOrder][ScaleBinSize
  // + 1][AngBinSize][ExtMomBinSize];

  // initialize interaction table
  for (int inin = 0; inin < AngBinSize; ++inin)
    for (int qIndex = 0; qIndex < ExtMomBinSize; ++qIndex) {
      for (int tIndex = 0; tIndex < TauBinSize; ++tIndex) {
        double k = Index2Mom(qIndex);
        double t = Index2Tau(tIndex);
        // EffInter(inin, qIndex, tIndex) = 8.0 * PI / (k * k + Para.Mass2) /
        //                                  (1 + pow(t * 10.0, 2) * 10.0 / PI);
        EffInter(inin, qIndex, tIndex) = 0.0;
        for (int order = 0; order < MaxOrder; ++order) {
          if (order == 0)
            DiffInter(order, inin, qIndex, tIndex) =
                EffInter(inin, qIndex, tIndex);
          else
            DiffInter(order, inin, qIndex, tIndex) = 0.0;
        }
      }
    }
}

double &verQTheta::EffInter(int Angle, int ExtQ, int Tau) {
  return EffInteraction[Angle * AngleIndex + ExtQ * QIndex + Tau];
}

double &verQTheta::DiffInter(int Order, int Angle, int ExtQ, int Tau) {
  return DiffInteraction[Order * OrderIndex + Angle * AngleIndex +
                         ExtQ * QIndex + Tau];
}

double verQTheta::Interaction(const momentum &InL, const momentum &InR,
                              const momentum &Transfer, double Tau,
                              int VerType) {

  double k = Transfer.norm();
  if (VerType == 0) {
    // return 8.0 * PI / (k * k + Para.Mass2) /
    //        (1 + pow(Tau * 10.0, 2) * 10.0 / PI);
    return -8.0 * PI / (k * k + Para.Mass2) / Para.Beta;
    // return 8.0 * PI / (k * k + Para.Mass2) *
    //        (0.5 / (1 + pow(abs(Tau) * 10.0, 2) * 10.0 / PI) +
    //         0.5 / (1 + pow((Para.Beta - abs(Tau)) * 10.0, 2) * 10.0 / PI));
    // return 1.0 / Para.Beta;
  } else if (VerType == 1) {
    if (k < Para.MaxExtMom) {
      if (Tau < 0.0)
        Tau += Para.Beta;
      int AngleIndex = Angle2Index(Angle2D(InL, InR), AngBinSize);
      int TauIndex = Tau2Index(Tau);
      return EffInter(AngleIndex, Mom2Index(k), TauIndex);
      // double Upper = EffInter(AngleIndex, Mom2Index(k), Tau2Index(Tau));
      // double Lower = EffInter(AngleIndex, Mom2Index(k), Tau2Index(Tau));
      // double UpperTau = Index2Tau(TauIndex + 1);
      // double LowerTau = Index2Tau(TauIndex);
      // // cout << Upper << " : " << Lower << endl;
      // return Lower + (Upper - Lower) / (UpperTau - LowerTau) * (Tau -
      // LowerTau);
    } else {
      return 0.0;
    }
  } else if (VerType == -1) {
    return 1.0;
  } else if (VerType == -2) {
    // return exp(-Transfer.norm() / Para.Kf);
    return 1.0;
    // return 8.0 * PI / (k * k + Para.Mass2) *
    //        (0.5 / (1 + pow(abs(Tau) * 10.0, 2) * 10.0 / PI) +
    //         0.5 / (1 + pow((Para.Beta - abs(Tau)) * 10.0, 2) * 10.0 /
    //         PI));
    // return 0.5 / (1 + pow(abs(Tau) * 10.0, 2) * 10.0 / PI) +
    //        0.5 / (1 + pow((Para.Beta - abs(Tau)) * 10.0, 2) * 10.0 / PI);
  } else {
    ABORT("VerType can not be " << VerType);
  }
}

void verQTheta::Measure(const momentum &InL, const momentum &InR,
                        const int QIndex, int Order, double *Tau,
                        double *Weight, double WeightFactor) {
  // cout << Order << ", " << DiagNum << endl;
  if (Order == 0) {
    Normalization += Weight[0] * WeightFactor;
    // Normalization += WeightFactor;
  } else {
    double Factor = 1.0 / pow(2.0 * PI, 2 * Order);
    int AngleIndex = Angle2Index(Angle2D(InL, InR), AngBinSize);
    for (int tIndex = 0; tIndex < (Order + 1) * 2; ++tIndex) {
      double dTau = Tau[tIndex] - Tau[0];
      if (dTau < 0.0)
        dTau += Para.Beta;
      // } else {
      int tBin = Tau2Index(dTau);
      // cout << AngleIndex << endl;
      // cout << InL[0] << "," << InL[1] << endl;
      // cout << InR[0] << "," << InR[1] << endl;
      // cout << "angle: " << Angle2D(InL, InR) << endl;
      DiffInter(Order, AngleIndex, QIndex, tBin) +=
          Weight[tIndex] * WeightFactor / Para.dAngleTable[AngleIndex] /
          (Para.Beta / TauBinSize) * Factor;

      // DiffInter(Order, AngleIndex, QIndex, tIndex) +=
      //     WeightFactor / Para.dAngleTable[AngleIndex] /
      //     (Para.Beta / TauBinSize) * Factor;
      // }
      DiffInter(0, AngleIndex, QIndex, tBin) +=
          Weight[tIndex] * WeightFactor / Para.dAngleTable[AngleIndex] * Factor;
    }
  }
  return;
}

void verQTheta::Update(double Ratio) {
  for (int angle = 0; angle < AngBinSize; ++angle)
    for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
      for (int tindex = 0; tindex < TauBinSize; ++tindex) {
        double OldValue = EffInter(angle, qindex, tindex);
        double NewValue = 0.0;
        for (int order = 1; order < MaxOrder; ++order) {
          // for (int order = 1; order < 2; ++order) {
          NewValue += DiffInter(order, angle, qindex, tindex) / Normalization *
                      PhyWeight;
        }
        EffInter(angle, qindex, tindex) =
            OldValue * (1 - Ratio) + NewValue * Ratio;
      }
}

void verQTheta::Save() {

  for (int order = 0; order < 4; order++) {
    string FileName = fmt::format("vertex{0}_pid{1}.dat", order, Para.PID);
    ofstream VerFile;
    VerFile.open(FileName, ios::out | ios::trunc);
    if (VerFile.is_open()) {
      VerFile << fmt::sprintf(
          "#PID:%d, Type:%d, rs:%.3f, Beta: %.3f, Step: %d\n", Para.PID,
          Para.ObsType, Para.Rs, Para.Beta, Para.Counter);
      VerFile << "# TauTable: ";
      for (int tau = 0; tau < TauBinSize; ++tau)
        VerFile << Index2Tau(tau) << " ";
      VerFile << endl;
      VerFile << "# AngleTable: ";
      for (int angle = 0; angle < AngBinSize; ++angle)
        VerFile << Para.AngleTable[angle] << " ";
      VerFile << endl;
      VerFile << "# ExtMomBinTable: ";
      for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
        VerFile << Para.ExtMomTable[qindex][0] << " ";
      VerFile << endl;

      for (int angle = 0; angle < AngBinSize; ++angle)
        for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
          for (int tindex = 0; tindex < TauBinSize; ++tindex)
            VerFile << DiffInter(order, angle, qindex, tindex) / Normalization *
                           PhyWeight
                    << "  ";
      VerFile.close();
    } else {
      LOG_WARNING("Polarization for PID " << Para.PID << " fails to save!");
    }
  }

  string FileName = fmt::format("vertex_pid{0}.dat", Para.PID);
  ofstream VerFile;
  VerFile.open(FileName, ios::out | ios::trunc);
  if (VerFile.is_open()) {
    VerFile << fmt::sprintf("#PID:%d, Type:%d, rs:%.3f, Beta: %.3f, Step: %d\n",
                            Para.PID, Para.ObsType, Para.Rs, Para.Beta,
                            Para.Counter);
    VerFile << "# TauTable: ";
    for (int tau = 0; tau < TauBinSize; ++tau)
      VerFile << Index2Tau(tau) << " ";
    VerFile << endl;
    VerFile << "# AngleTable: ";
    for (int angle = 0; angle < AngBinSize; ++angle)
      VerFile << Para.AngleTable[angle] << " ";
    VerFile << endl;
    VerFile << "# ExtMomBinTable: ";
    for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
      VerFile << Para.ExtMomTable[qindex][0] << " ";
    VerFile << endl;

    for (int angle = 0; angle < AngBinSize; ++angle)
      for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
        for (int tindex = 0; tindex < TauBinSize; ++tindex)
          VerFile << EffInter(angle, qindex, tindex) << "  ";
    VerFile.close();
  } else {
    LOG_WARNING("Polarization for PID " << Para.PID << " fails to save!");
  }
}

void verQTheta::ClearStatis() {
  Normalization = 1.0e-10;
  for (int inin = 0; inin < AngBinSize; ++inin)
    for (int qIndex = 0; qIndex < ExtMomBinSize; ++qIndex) {
      double k = Index2Mom(qIndex);
      for (int order = 0; order < MaxOrder; ++order)
        for (int tindex = 0; tindex < TauBinSize; ++tindex)
          DiffInter(order, inin, qIndex, tindex) = 0.0;
    }
}

fermi::fermi() {
  UpperBound = 5.0 * Para.Ef;
  LowerBound = 0.0;
  DeltaK = UpperBound / MAXSIGMABIN;
  UpperBound2 = 1.2 * Para.Ef;
  LowerBound2 = 0.8 * Para.Ef;
  DeltaK2 = UpperBound2 / MAXSIGMABIN;
  if (Para.SelfEnergyType == FOCK)
    BuildFockSigma();
}

double fermi::Fock(double k) {
  // warning: this function only works for T=0!!!!
  double l = sqrt(Para.Mass2);
  double kF = Para.Kf;
  double fock = 1.0 + l / kF * atan((k - kF) / l);
  fock -= l / kF * atan((k + kF) / l);
  fock -= (l * l - k * k + kF * kF) / 4.0 / k / kF *
          log((l * l + (k - kF) * (k - kF)) / (l * l + (k + kF) * (k + kF)));
  fock *= (-2.0 * kF) / PI;

  double shift = 1.0 - l / kF * atan(2.0 * kF / l);
  shift -= l * l / 4.0 / kF / kF * log(l * l / (l * l + 4.0 * kF * kF));
  shift *= (-2.0 * kF) / PI;

  return fock - shift;
}

double fermi::BuildFockSigma() {
  ASSERT_ALLWAYS(D == 3, "The Fock self energy is for 3D!");
  double fock, k;
  for (int i = 0; i < MAXSIGMABIN; ++i) {
    // k: (0^+, UpBound^-)
    // i=0 ==> k==0.5*DeltaK
    // i=MAXSIGMABIN-1 ==> k==(MAXSIGMABIN-0.5)*DeltaK
    k = (i + 0.5) * DeltaK + LowerBound;
    Sigma[i] = Fock(k);
    if (i > 0 && k <= LowerBound2 && k >= UpperBound2) {
      ASSERT_ALLWAYS(
          Equal(Sigma[i - 1], Sigma[i], 5.0e-5),
          fmt::format("Fock are not accurate enough! At k={0}: {1} vs {2}\n", k,
                      Sigma[i - 1], Sigma[i]));
    }
    // cout << k << " : " << Sigma[i] << " vs " << Fock(k) << endl;
  }

  for (int i = 0; i < MAXSIGMABIN; ++i) {
    // k: (0^+, UpBound^-)
    // i=0 ==> k==0.5*DeltaK
    // i=MAXSIGMABIN-1 ==> k==(MAXSIGMABIN-0.5)*DeltaK
    k = (i + 0.5) * DeltaK2 + LowerBound2;
    Sigma2[i] = Fock(k);
    if (i > 0) {
      ASSERT_ALLWAYS(Equal(Sigma2[i - 1], Sigma2[i], 5.0e-5),
                     fmt::format("The 2rd level Fock are not accurate enough!"
                                 "level! At k={0}: {1} vs {2}\n",
                                 k, Sigma2[i - 1], Sigma2[i]));
    }
    // cout << k << " : " << Sigma[i] << " vs " << Fock(k) << endl;
  }
};

double fermi::FockSigma(const momentum &Mom) {
  double k = Mom.norm(); // bare propagator
  double fock;
  if (k >= LowerBound2 && k < UpperBound2) {
    int i = (k - LowerBound2) / DeltaK2;
    fock = Sigma2[i];
  } else if ((k >= LowerBound && k < LowerBound2) ||
             (k >= UpperBound2 && k < UpperBound)) {
    int i = (k - LowerBound) / DeltaK;
    fock = Sigma[i];
  } else {
    fock = Fock(k);
  }
  // ASSERT_ALLWAYS(
  //     Equal(fock, Fock(k), 5.0e-5),
  //     fmt::format("Fock are not accurate enough! At k={0}: {1} vs {2}\n",
  //     k,
  //                 fock, Fock(k)));
  return fock + k * k;
}

double fermi::PhyGreen(double Tau, const momentum &Mom, int GType,
                       double Scale) {
  // return 1.0;
  // if tau is exactly zero, set tau=0^-
  double green, Ek, kk, k;
  if (Tau == 0.0) {
    return EPS;
  }
  // equal time green's function
  if (GType == 1)
    Tau = -1.0e-10;

  double s = 1.0;
  if (Tau < 0.0) {
    Tau += Para.Beta;
    s = -s;
  } else if (Tau >= Para.Beta) {
    Tau -= Para.Beta;
    s = -s;
  }

  k = Mom.norm();
  if (Para.SelfEnergyType == BARE)
    Ek = k * k; // bare propagator
  else if (Para.SelfEnergyType == FOCK)
    Ek = FockSigma(Mom); // Fock diagram dressed propagator
  else
    ABORT("Green function is not implemented!");

  //// enforce an UV cutoff for the Green's function ////////
  // if(Ek>8.0*EF) then
  //   PhyGreen=0.0
  //   return
  // endif

  double x = Para.Beta * (Ek - Para.Mu) / 2.0;
  double y = 2.0 * Tau / Para.Beta - 1.0;
  if (x > 100.0)
    green = exp(-x * (y + 1.0));
  else if (x < -100.0)
    green = exp(x * (1.0 - y));
  else
    green = exp(-x * y) / (2.0 * cosh(x));

  green *= s;

  // if (std::isnan(green))
  //   ABORT("Step:" << Para.Counter << ", Green is too large! Tau=" << Tau
  //                 << ", Ek=" << Ek << ", Green=" << green << ", Mom"
  //                 << ToString(Mom));
  return green;
}

double fermi::Green(double Tau, const momentum &Mom, spin Spin, int GType,
                    double Scale) {
  double green;
  if (GType >= 0) {
    green = PhyGreen(Tau, Mom, GType, Scale);
  } else if (GType == -1) {
    // green = PhyGreen(Tau, Mom, Scale);
    green = 1.0;

  } else if (GType == -2) {
    // Lower Scale Green's function
    Scale -= 1;
    green = PhyGreen(Tau, Mom, GType, Scale);
  } else {
    ABORT("GType " << GType << " has not yet been implemented!");
    // return FakeGreen(Tau, Mom);
  }
  return green;
}

verfunc::verfunc() {
  // test angle utility

  // TODO: implement D=3
  if (D == 3)
    return;

  // _TestAngle2D();
  // _TestAngleIndex();

  // // initialize UV ver4 table
  // momentum KInL = {1.0, 0.0};
  // for (int inin = 0; inin < InInAngBinSize; ++inin)
  //   for (int inout = 0; inout < InOutAngBinSize; ++inout) {
  //     double AngleInIn = Index2Angle(inin, InInAngBinSize);
  //     double AngleInOut = Index2Angle(inout, InOutAngBinSize);
  //     momentum KInR = {cos(AngleInIn), sin(AngleInIn)};
  //     momentum KOutL = {cos(AngleInOut), sin(AngleInOut)};
  //     momentum KOutR = KInL + KInR - KOutL;
  //     Ver4AtUV[inin][inout] = (KInL - KInR).dot(KOutL - KOutR);
  //   }
}

void verfunc::Vertex4(const momentum &InL, const momentum &InR,
                      const momentum &OutL, const momentum &OutR,
                      int Ver4TypeDirect, int Ver4TypeExchange, double &Direct,
                      double &Exchange) {
  if (Ver4TypeDirect != 0 || Ver4TypeExchange != 0)
    ABORT("Ver4Type is only implemented for 0!");

  /**************   Yokawar Interaction ************************/
  Direct = 8.0 * PI / ((OutL - InL).squaredNorm() + Para.Mass2);
  Exchange = 8.0 * PI / ((OutR - InL).squaredNorm() + Para.Mass2);

  /**************   Generic Interaction ************************/
}

double diag::Index2Mom(const int &Index) {
  return (Index + 0.5) / ExtMomBinSize * Para.MaxExtMom;
};

int diag::Mom2Index(const double &K) {
  return int(K / Para.MaxExtMom * ExtMomBinSize);
};

double diag::Angle2D(const momentum &K1, const momentum &K2) {
  // Returns the angle in radians between vectors 'K1' and 'K2'
  // double dotp = K1.dot(K2);
  double dotp = K1[0] * K2[0] + K1[1] * K2[1];
  double det = K1[0] * K2[1] - K1[1] * K2[0];
  double Angle2D = atan2(det, dotp);

  // cout<<endl;
  // cout << K1[0] << "," << K1[1] << endl;
  // cout << K2[0] << "," << K2[1] << endl;
  // cout << "dotp:" << dotp << endl;
  // cout << "det:" << det << endl;
  // cout << "angle:" << Angle2D << endl;
  if (Angle2D < 0)
    Angle2D += 2.0 * PI;
  // cout << "angleadjusted:" << Angle2D << endl;
  // cout << endl;
  return Angle2D;
}

double diag::Index2Angle(const int &Index, const int &AngleNum) {
  // Map index [0...AngleNum-1] to the theta range [0.0, 2*pi)
  return Index * 2.0 * PI / AngleNum;
}

int diag::Angle2Index(const double &Angle, const int &AngleNum) {
  // Map theta range  [0.0, 2*pi) to index [0...AngleNum-1]
  double dAngle = 2.0 * PI / AngleNum;
  if (Angle >= 2.0 * PI - dAngle / 2.0 || Angle < dAngle / 2.0)
    return 0;
  else
    return int(Angle / dAngle + 0.5);
}

double diag::Index2Scale(const int &Index) {
  return Index * Para.UVScale / ScaleBinSize;
}
int diag::Scale2Index(const double &Scale) {
  return int((Scale / Para.UVScale) * ScaleBinSize);
}

int diag::Tau2Index(const double &Tau) {
  return int((Tau / Para.Beta) * TauBinSize);
}

double diag::Index2Tau(const int &Index) {
  return (Index + 0.5) * Para.Beta / TauBinSize;
}

void diag::_TestAngle2D() {
  // Test Angle functions
  momentum K1 = {1.0, 0.0};
  momentum K2 = {1.0, 0.0};

  ASSERT_ALLWAYS(
      abs(Angle2D(K1, K2)) < 1.e-7,
      fmt::format("Angle between K1 and K2 are not zero! It is {:.13f}",
                  Angle2D(K1, K2)));

  K1 = {1.0, 0.0};
  K2 = {-1.0, 0.0};
  ASSERT_ALLWAYS(
      abs(Angle2D(K1, K2) - PI) < 1.e-7,
      fmt::format("Angle between K1 and K2 are not Pi! Instead, it is {:.13f}",
                  Angle2D(K1, K2)));

  K1 = {1.0, 0.0};
  K2 = {1.0, -EPS};
  ASSERT_ALLWAYS(
      abs(Angle2D(K1, K2) - 2.0 * PI) < 1.e-7,
      fmt::format("Angle between K1 and K2 are not 2.0*Pi! It is {:.13f}",
                  Angle2D(K1, K2)));
}

void diag::_TestAngleIndex() {
  // Test Angle functions
  int AngleNum = 64;
  ASSERT_ALLWAYS(abs(Index2Angle(0, AngleNum) - 0.0) < 1.0e-10,
                 "Angle for index 0 should be zero!");

  ASSERT_ALLWAYS(abs(Index2Angle(AngleNum - 1, AngleNum) -
                     (2.0 * PI * (1.0 - 1.0 / AngleNum))) < 1.0e-10,
                 "Angle for index AngleNum should be 2.0*pi-0^+!");

  ASSERT_ALLWAYS(Angle2Index(0.0, AngleNum) == 0,
                 "Angle zero should have index 0!");
  ASSERT_ALLWAYS(
      Angle2Index(2.0 * PI * (1.0 - 0.5 / AngleNum) + EPS, AngleNum) == 0,
      "Angle 2*pi-pi/AngleNum should have index 1!");

  ASSERT_ALLWAYS(Angle2Index(2.0 * PI * (1.0 - 0.5 / AngleNum) - EPS,
                             AngleNum) == AngleNum - 1,
                 "Angle 2*pi-pi/AngleNum-0^+ should have index AngleNum!");
}
