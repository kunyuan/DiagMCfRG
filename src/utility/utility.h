//---------------------------------------------------------------------------
// �Nikolai V. Shokhirev, 2004-2008  <nikolai@shokhirev.com>
// http://www.shokhirev.com/nikolai.html
// Reduced demo version
//---------------------------------------------------------------------------

#ifndef MathUtils_H
#define MathUtils_H

#include "abort.h"
#include "stdio.h"
#include <math.h>
#include <sstream>
#include <string>

//---------------------------------------------------------------------------
const double eps0 = 1.0e-9;

// FORTRAN abs
double fabs(double x); // { return ( (x >= 0.0) ? x : -x); }

// FORTRAN iabs
int iabs(int x); // { return ( (x >= 0.0) ? x : -x); }
// more functions ...

// float iszero
bool Zero(double x, double eps = eps0);

// float equal
bool Equal(double x1, double x2, double eps = eps0);
bool Equal(uint x1, uint x2, double eps = eps0);
bool Equal(int x1, int x2, double eps = eps0);
template <typename T>
bool Equal(const T *x1, const T *x2, uint num, double eps = eps0) {
  for (uint i = 0; i < num; i++)
    if (!Equal(x1[i], x2[i], eps0))
      return false;
  return true;
}

template <typename T> std::string ToString(const T &value) {
  std::ostringstream oss;
  oss << value;
  return oss.str();
}

template <typename T> std::string ToString(const T *array, size_t Num) {
  std::ostringstream oss;
  oss << "[";
  for (size_t i = 0; i < Num; i++)
    oss << " " << array[i] << ",";
  oss << "]";
  return oss.str();
}

/* Convert double to string with specified number of places after the decimal
   and left padding. */
std::string ToString(const double x, const int width, const int decDigits = 5);
/*! Center-aligns string within a field of width w. Pads with blank spaces
    to enforce alignment. */
std::string Center(const std::string s, const int w);
/* Right-aligns string within a field of width w. Pads with blank spaces
 to enforce alignment. */
std::string Right(const std::string s, const int w);
/*! Left-aligns string within a field of width w. Pads with blank spaces
 to enforce alignment. */
std::string Left(const std::string s, const int w);

template <typename T> void AssignFromTo(T *source, T *target, int size) {
  for (int i = 0; i < size; i++)
    target[i] = source[i];
}

template <typename T> void InitialArray(T *target, T t, const int &size) {
  for (int i = 0; i < size; i++)
    target[i] = t;
}

bool CleanFile(const std::string &FileName);

bool DoesFileExist(const std::string &FileName);

#define CHECKNULL(source)                                                      \
  {                                                                            \
    if ((source) == nullptr)                                                   \
      ABORT(#source << " is nullptr!");                                        \
  }

std::string ProgressBar(double progress);

// template <typename... TS>
// static std::string Format(std::string format, TS... args) {
//   char format_buffer[512];
//   int length = sprintf(format_buffer, format.c_str(), args...);
//   ASSERT_ALLWAYS(length > 0, "Failed to construct the string!");
//   return std::string(format_buffer, length);
// }

#endif