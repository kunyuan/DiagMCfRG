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

  PhyWeight = ExtMomBinSize * 2.0 * PI * Para.Kf * 4.0;

  // initialize interaction table
  for (int scale = 0; scale < ScaleBinSize + 1; ++scale) {
    for (int inin = 0; inin < AngBinSize; ++inin)
      for (int qIndex = 0; qIndex < ExtMomBinSize; ++qIndex) {
        double k = Index2Mom(qIndex);
        EffInteraction[scale][inin][qIndex] = 8.0 * PI / (k * k + Para.Mass2);
        // EffInteraction[scale][inin][qIndex] = 0.0;
        for (int order = 0; order < MaxOrder; ++order)
          DiffInteraction[order][scale][inin][qIndex] = 0.0;
      }
  }
}

double verQTheta::Interaction(const momentum &InL, const momentum &InR,
                              const momentum &Transfer, int VerType,
                              double Scale) {
  if (VerType >= 0) {

    double k = Transfer.norm();
    // return 8.0 * PI / (k * k + Para.Mass2);
    // return 1.0;

    int AngleIndex = Angle2Index(Angle2D(InL, InR), AngBinSize);
    if (k < Para.MaxExtMom && Scale < Para.UVScale) {
      // return 8.0 * PI /
      // (k * k + Para.Mass2 +
      // 0.159 * (1 - sqrt(1 - 4.0 * Para.Kf * Para.Kf / k / k)));
      int ScaleIndex = Scale2Index(Scale);
      double Upper = EffInteraction[ScaleIndex + 1][AngleIndex][Mom2Index(k)];
      double Lower = EffInteraction[ScaleIndex][AngleIndex][Mom2Index(k)];
      double UpperScale = Para.ScaleTable[ScaleIndex + 1];
      double LowerScale = Para.ScaleTable[ScaleIndex];
      return Lower +
             (Upper - Lower) / (UpperScale - LowerScale) * (Scale - LowerScale);
    } else
      return 8.0 * PI / (k * k + Para.Mass2);

  } else if (VerType == -1) {
    return 1.0;
  } else if (VerType == -2) {
    // return exp(-Transfer.norm() / Para.Kf);
    return exp(-Scale / Para.Kf / 4.0);
  } else {
    ABORT("VerType can not be " << VerType);
  }
}

void verQTheta::Measure(const momentum &InL, const momentum &InR,
                        const int QIndex, double Scale, int Order,
                        double WeightFactor) {
  if (Order == 0) {
    Normalization += WeightFactor;
  } else {
    int AngleIndex = Angle2Index(Angle2D(InL, InR), AngBinSize);
    int ScaleIndex = Scale2Index(Scale);
    // cout << AngleIndex << endl;
    // cout << InL[0] << "," << InL[1] << endl;
    // cout << InR[0] << "," << InR[1] << endl;
    // cout << "angle: " << Angle2D(InL, InR) << endl;
    DiffInteraction[Order][ScaleIndex][AngleIndex][QIndex] +=
        WeightFactor / Para.dAngleTable[AngleIndex];
  }
  return;
}

void verQTheta::Update(double Ratio) {
  for (int scale = ScaleBinSize - 1; scale >= 0; --scale)
    for (int angle = 0; angle < AngBinSize; ++angle)
      for (int qindex = 0; qindex < ExtMomBinSize; ++qindex) {
        double OldValue = EffInteraction[scale][angle][qindex];
        double NewValue = EffInteraction[scale + 1][angle][qindex];
        // double NewValue = 0.0;
        for (int order = 0; order < MaxOrder; ++order) {
          NewValue += DiffInteraction[order][scale][angle][qindex] /
                      Normalization * PhyWeight / 40.0;
        }
        EffInteraction[scale][angle][qindex] =
            OldValue * (1 - Ratio) + NewValue * Ratio;
      }
}

void verQTheta::Save() {
  string FileName = fmt::format("vertex_pid{0}.dat", Para.PID);
  ofstream VerFile;
  VerFile.open(FileName, ios::out | ios::trunc);
  if (VerFile.is_open()) {
    VerFile << fmt::sprintf("#PID:%d, Type:%d, rs:%.3f, Beta: %.3f, Step: %d\n",
                            Para.PID, Para.ObsType, Para.Rs, Para.Beta,
                            Para.Counter);
    VerFile << "# ScaleTable: ";
    for (int scale = 0; scale < ScaleBinSize + 1; ++scale)
      VerFile << Para.ScaleTable[scale] << " ";
    VerFile << endl;
    VerFile << "# AngleTable: ";
    for (int angle = 0; angle < AngBinSize; ++angle)
      VerFile << Para.AngleTable[angle] << " ";
    VerFile << endl;
    VerFile << "# ExtMomBinTable: ";
    for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
      VerFile << Para.ExtMomTable[qindex][0] << " ";
    VerFile << endl;

    for (int scale = 0; scale < ScaleBinSize + 1; ++scale)
      for (int angle = 0; angle < AngBinSize; ++angle)
        for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
          VerFile << EffInteraction[scale][angle][qindex] << "  ";
    VerFile.close();
  } else {
    LOG_WARNING("Polarization for PID " << Para.PID << " fails to save!");
  }
}

void verQTheta::ClearStatis() {
  Normalization = 1.0e-10;
  for (int scale = 0; scale < ScaleBinSize + 1; ++scale) {
    for (int inin = 0; inin < AngBinSize; ++inin)
      for (int qIndex = 0; qIndex < ExtMomBinSize; ++qIndex) {
        double k = Index2Mom(qIndex);
        for (int order = 0; order < MaxOrder; ++order)
          DiffInteraction[order][scale][inin][qIndex] = 0.0;
      }
  }
}

void verQTheta::ResetIRScale(int IRScaleBin) {
  for (int scale = 0; scale < IRScaleBin; ++scale) {
    for (int inin = 0; inin < AngBinSize; ++inin)
      for (int qIndex = 0; qIndex < ExtMomBinSize; ++qIndex) {
        double k = Index2Mom(qIndex);
        EffInteraction[scale][inin][qIndex] =
            EffInteraction[IRScaleBin][inin][qIndex];
      }
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

cmplx fermi::GreenFreq(double Freq, const momentum &Mom, int GType,
                       double Scale) {
  cmplx green;
  double k = Mom.norm();
  green = -1.0 / (1i * Freq - k * k);
  if (Para.Type == RG) {
    double kScale = Scale;
    double dK2 = (k - Para.Kf) * (k - Para.Kf);

    if (GType == 2)
      green *= -dK2 * 4.0 * pow(kScale, 3) / pow((dK2 + pow(kScale, 4)), 2);
    else
      green *= dK2 / (dK2 + pow(kScale, 4));

    // if (GType == 2)
    //   green *= -dK2 / (dK2 + kScale) / (dK2 + kScale);
    // else
    //   green *= dK2 / (dK2 + kScale);
  }
  return green;
}

double fermi::PhyGreen(double Tau, const momentum &Mom, int GType,
                       double Scale) {
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

  // if (Debug) { //   cout << "Tau=" << Tau << endl;
  //   cout << "Counter" << Para.Counter << endl;
  //   cout << "Tau=" << Tau << endl;
  //   cout << "Mom=(" << Mom[0] << ", " << Mom[1] << ")" << endl;
  //   cout << "Mom2=" << Mom[0] * Mom[0] + Mom[1] * Mom[1] << endl;
  //   cout << "Ek=" << Ek << endl;
  //   cout << "x=" << x << endl;
  //   cout << "y=" << y << endl;
  //   cout << "Green=" << green << endl << endl;
  // }

  // cout << "x: " << x << ", y: " << y << ", G: " << green << endl;
  // cout << "G: " << green << endl;

  if (Para.Type == RG) {
    // double kScale = Para.ScaleTable[Scale];
    double kScale = Scale;
    // double dK2 = (k - Para.Kf) * (k - Para.Kf) / kScale / kScale;
    // double expFactor = 0.0;
    // if (dK2 < 100.0)
    //   expFactor = exp(-dK2);
    // green *= (1 - expFactor);
    // if (GType == 2)
    //   green *= -expFactor * 2.0 * dK2 / kScale;
    double dK2 = (k - Para.Kf) * (k - Para.Kf);
    // if (GType == 2)
    //   green *= -dK2 * 2.0 * kScale / (dK2 + kScale * kScale) /
    //            (dK2 + kScale * kScale);
    // else
    //   green *= dK2 / (dK2 + kScale * kScale);

    if (GType == 2)
      green *= -dK2 * 4.0 * pow(kScale, 3) / pow((dK2 + pow(kScale, 4)), 2);
    else
      green *= dK2 / (dK2 + pow(kScale, 4));

    // if (GType == 2)
    //   green *= -dK2 / (dK2 + kScale) / (dK2 + kScale);
    // else
    //   green *= dK2 / (dK2 + kScale);
  }

  if (std::isnan(green))
    ABORT("Step:" << Para.Counter << ", Green is too large! Tau=" << Tau
                  << ", Ek=" << Ek << ", Green=" << green << ", Mom"
                  << ToString(Mom));
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
