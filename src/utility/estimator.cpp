//
//  Estimate.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/14/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "estimator.h"
#include "abort.h"
#include "global.h"
#include <algorithm>
#include <iostream>

using namespace std;

template <typename T> EstimateClass<T>::EstimateClass() : Mean(), Error() {}

/**
 *  Estimate constructor of type T
 *
 *  @param mean __reference__ of a T number as the mean value
 *  @param error __reference__ of a T number as the error
 *
 */
template <typename T>
EstimateClass<T>::EstimateClass(const T &mean, const T &error)
    : Mean(mean), Error(error) {}
template <> double EstimateClass<double>::RelativeError() {
  return Error / Mean;
}

/**
 *  \brief Pretty output of Complex Estimate
 */

ostream &operator<<(ostream &os, const EstimateClass<double> &e) {
  os.setf(ios::showpoint);
  os << e.Mean << "+/-" << e.Error;
  return os;
}

template class EstimateClass<double>;

/**
 *  Normalization factor is initialized as 1.0 when Estimator is constructed.
 *
 */
template <typename T> Estimator<T>::Estimator() {
  Name = "Temp";
  ClearStatistics();
}

template <typename T> Estimator<T>::Estimator(string name) {
  Name = name;
  ClearStatistics();
}

template <typename T> void Estimator<T>::ClearStatistics() {
  _interval = 1;
  _counter = 0;
  _end = 0;
  _accumulator = 0.0;
  _ratio = 1.0;
  _norm = 1.0;
  // small enough non-zero number to avoid NAN
}

template <typename T> void Estimator<T>::SqueezeStatistics(double factor) {
  ASSERT_ALLWAYS(factor > 0.0, "factor=" << factor << "<=0!");
  _accumulator /= factor;
  _norm /= factor;
}
/**
 *  \brief Using statistics from ThrowRatio*100% to 100% to estimate the error
 * bar
 */
const double ThrowRatio = 1.0 / 3;
template <> void Estimator<double>::_update() {
  _value.Mean = _accumulator / _norm;
  uint size = _end;
  if (size == 0)
    return;
  double Min = MAXREAL, Max = MINREAL;
  int MinIndex = 0, MaxIndex = 0;
  for (uint i = size * ThrowRatio; i < size; i++) {
    if (Min > _history[i]) {
      Min = _history[i];
      MinIndex = i;
    }
    if (Max < _history[i]) {
      Max = _history[i];
      MaxIndex = i;
    }
  }
  _value.Error = fabs(Max - Min) / 2.0;
  _ratio = (MaxIndex - MinIndex) / (double)size * (1.0 - ThrowRatio);
}

// template <>
// void Estimator<Complex>::_update()
// {
//     _value.Mean = _accumulator / _norm;

//     uint size = _end;
//     if (size == 0)
//         return;
//     Complex Min(MaxReal, MaxReal), Max(MinReal, MinReal);
//     int MinIndexRe = 0, MaxIndexRe = 0;
//     int MinIndexIm = 0, MaxIndexIm = 0;
//     for (uint i = size * ThrowRatio; i < size; i++) {
//         if (Min.Re > _history[i].Re) {
//             Min.Re = _history[i].Re;
//             MinIndexRe = i;
//         }
//         if (Max.Re < _history[i].Re) {
//             Max.Re = _history[i].Re;
//             MaxIndexRe = i;
//         }
//         if (Min.Im > _history[i].Im) {
//             Min.Im = _history[i].Im;
//             MinIndexIm = i;
//         }
//         if (Max.Im < _history[i].Im) {
//             Max.Im = _history[i].Im;
//             MaxIndexIm = i;
//         }
//     }
//     _value.Error.Re = fabs(Max.Re - Min.Re) / 2.0;
//     _value.Error.Im = fabs(Max.Im - Min.Im) / 2.0;
//     if (MaxIndexRe - MinIndexRe < MaxIndexIm - MinIndexIm)
//         _ratio = (MaxIndexIm - MinIndexIm) / (real)size * (1.0 - ThrowRatio);
//     else
//         _ratio = (MaxIndexRe - MinIndexRe) / (real)size * (1.0 - ThrowRatio);
// }

template <typename T> void Estimator<T>::Measure(const T &t) {
  _accumulator += t;
  _norm += 1.0;
}

template <typename T> void Estimator<T>::AddStatistics() {
  if (_end == SIZE) {
    for (uint i = 0; i < _end / 2; i++)
      _history[i] = _history[2 * i];
    _end /= 2;
    _interval *=
        2; // perform measurement for every _interval calls of AddStatistics
  } else if (_end < SIZE) {
    _counter++;
    if (_counter >= _interval) {
      _history[_end] = _accumulator / _norm;
      _end++;
      _counter = 0;
    }
  } else
    ABORT("Exceed history memory!");
}

template <typename T> EstimateClass<T> Estimator<T>::Estimate() {
  _update();
  return _value;
}

template <typename T> T Estimator<T>::Value() { return _accumulator; }

template <typename T> double Estimator<T>::Norm() { return _norm; }

template <typename T> double Estimator<T>::Ratio() { return _ratio; }

template class Estimator<double>;

// bundle of estimators
template <typename T> void EstimatorBundle<T>::AddEstimator(string name) {
  _MakeSureKeyNotExists(name);
  _EstimatorVector.push_back(EstimatorT(name));
  _EstimatorMap[name] = _EstimatorVector.size() - 1;
}

/**
 *  \brief this function will give you a new copy of Estimator<T>, including a
 * __new__ Estimator<T>._history
 */
template <typename T>
void EstimatorBundle<T>::AddEstimator(const Estimator<T> &est) {
  _MakeSureKeyNotExists(est.Name);
  _EstimatorVector.push_back(est);
  _EstimatorMap[est.Name] = _EstimatorVector.size() - 1;
}

template <typename T>
bool EstimatorBundle<T>::_MakeSureKeyNotExists(string key) {
  if (_EstimatorMap.find(key) != _EstimatorMap.end()) {
    ABORT(key << " has already existed!");
    return true;
  }
  return false;
}

template <typename T> void EstimatorBundle<T>::AddStatistics() {
  for (int i = 0; i < HowMany(); i++)
    _EstimatorVector[i].AddStatistics();
}

template <typename T> int EstimatorBundle<T>::HowMany() {
  return (int)_EstimatorVector.size();
}

template <typename T> Estimator<T> &EstimatorBundle<T>::operator[](int index) {
  return _EstimatorVector[index];
}

template <typename T>
Estimator<T> &EstimatorBundle<T>::operator[](string name) {
  return _EstimatorVector[_EstimatorMap[name]];
}

/**
 *  \brief clear all statistics of the elements in the EstimatorBundle.
 *   __memory__ of the vector will not be freed!
 *
 */
template <typename T> void EstimatorBundle<T>::ClearStatistics() {
  for (auto &vector : _EstimatorVector)
    vector.ClearStatistics();
}

template <typename T>
void EstimatorBundle<T>::SqueezeStatistics(double factor) {
  for (auto &vector : _EstimatorVector)
    vector.SqueezeStatistics(factor);
}

template class EstimatorBundle<double>;