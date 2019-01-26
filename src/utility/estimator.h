//
//  Estimate.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/14/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __FeynCalc__Estimate__
#define __FeynCalc__Estimate__

#include <iosfwd>
#include <unordered_map>
#include <vector>

#define SIZE 100000

/**
 *  \brief estimate with mean value and standard error
 */
template <typename T> class EstimateClass {
public:
  EstimateClass();
  EstimateClass(const T &m, const T &e);
  T Mean;
  T Error;
  T RelativeError();
  friend std::ostream &operator<<(std::ostream &, const EstimateClass &);
};

/**
 *  \brief maintain a history of observables, where an Estimate object can be
 * calculated
 */

template <typename T> class Estimator {
private:
  uint _interval;
  uint _counter;
  uint _end;
  T _accumulator;
  double _norm;
  double _ratio;
  EstimateClass<T> _value;
  void _update();
  T _history[SIZE];

public:
  Estimator();
  Estimator(std::string);
  std::string Name;
  void Measure(const T &);
  void AddStatistics();
  T Value();
  double Norm();
  double Ratio();
  EstimateClass<T> Estimate();
  void ClearStatistics();
  void SqueezeStatistics(double factor);
};
/**
 *  \brief an extension to the std::vector<Estimator<T>> to contain
 * Estimator<T>
 */

template <typename T> class EstimatorBundle {
private:
  typedef Estimator<T> EstimatorT;
  std::vector<EstimatorT> _EstimatorVector;
  std::unordered_map<std::string, unsigned long> _EstimatorMap;
  bool _MakeSureKeyNotExists(std::string key);

public:
  void AddEstimator(const std::string);
  void AddEstimator(const EstimatorT &);
  void RemoveEstimator(const std::string);
  void AddStatistics();
  int HowMany();
  // if AllowFailure is set to true, then FromDict will not throw exception even
  // if some of the estimators fail to read
  EstimatorT &operator[](int);
  EstimatorT &operator[](std::string);
  void ClearStatistics();
  void SqueezeStatistics(double factor);
};

#endif /* defined(__Feynman_Simulator__Estimate__) */
