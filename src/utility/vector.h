
#ifndef __FeynCalc__vector__
#define __FeynCalc__vector__

#include <iosfwd>
#include <sstream>
#include <vector>
#include <math.h>
#include <initializer_list>
using namespace std;

template <typename T, int D>
class Vec
{
  private:
    T _Array[D];

  public:
    Vec()
    {
    }

    Vec(std::initializer_list<T> list)
    {
        int i = 0;
        for (auto p = list.begin(); p < list.end() && i < D; ++p)
        {
            _Array[i] = *p;
            i++;
        }
    }

    Vec(T value)
    {
        for (int j = 0; j < D; j++)
            _Array[j] = value;
    }

    Vec(T *value)
    {
        for (int j = 0; j < D; j++)
            _Array[j] = value[j];
    }

    void CopyToArray(T *target) const
    {
        for (int i = 0; i < D; i++)
            target[i] = _Array[i];
    }

    const T *begin() const
    {
        return _Array;
    }
    const T *end() const
    {
        return _Array + D;
    }
    uint size() const
    {
        return D;
    }

    T &operator[](int index)
    {
        return _Array[index];
    }

    const T &operator[](int index) const
    {
        return _Array[index];
    }

    Vec operator*(int i) const
    {
        Vec v;
        for (int j = 0; j < D; j++)
            v[j] = _Array[j] * i;
        return v;
    }

    Vec operator*(const double &i) const
    {
        Vec v;
        for (int j = 0; j < D; j++)
            v[j] = _Array[j] * i;
        return v;
    }

    Vec operator+(const Vec &v2) const
    {
        Vec v;
        for (int j = 0; j < D; j++)
            v[j] = _Array[j] + v2._Array[j];
        return v;
    }

    Vec operator-(const Vec &v2) const
    {
        Vec v;
        for (int j = 0; j < D; j++)
            v[j] = _Array[j] - v2._Array[j];
        return v;
    }

    Vec &operator+=(const Vec &v2)
    {
        for (int j = 0; j < D; j++)
            _Array[j] += v2._Array[j];
        return *this;
    }

    Vec &operator-=(const Vec &v2)
    {
        for (int j = 0; j < D; j++)
            _Array[j] -= v2._Array[j];
        return *this;
    }

    // string PrettyString(){
    //     std::ostringstream oss;
    //     oss << *this;
    //     return "(" + oss.str() + ")";
    // }

    // template <typename TT>
    // friend std::ostream& operator<<(std::ostream& os, const Vec<TT, D>&){
    //     for (int i = 0; i < D - 1; i++)
    //         os << v[i] << " ";
    //     os << v[D - 1];
    //     return os;
    // }

    // template <typename TT>
    // friend std::istream& operator>>(std::istream& is, Vec<TT, D>&){
    //     is >> v[0];
    //     char sep;
    //     for (int i = 1; i < D; i++) {
    //         is >> sep >> v[i];
    //         if (sep != " ")
    //             is.setstate(ios::failbit);
    //     }
    //     return is;
    // }

    // friend bool operator==(const Vec<int>&, const Vec<int>&);
    // friend bool operator==(const Vec<real>&, const Vec<real>&);
};
// bool operator!=(const Vec<int, D>& v1, const Vec<int, D>& v2);
// bool operator!=(const Vec<double, D>& v1, const Vec<double, D>& v2);

template <typename T, int D>
double dot(const Vec<T, D> &v1, const Vec<T, D> &v2)
{
    double sum = v1[0] * v2[0];
    for (int i = 1; i < D; i++)
    {
        sum += v1[i] * v2[i];
    }
    return sum;
}

template <typename T, int D>
double sum2(const Vec<T, D> &v1)
{
    double sum = v1[0] * v1[0];
    for (int i = 1; i < D; i++)
    {
        sum += v1[i] * v1[i];
    }
    return sum;
}

template <typename T, int D>
double norm2(const Vec<T, D> &v1)
{
    double sum = v1[0] * v1[0];
    for (int i = 1; i < D; i++)
    {
        sum += v1[i] * v1[i];
    }
    return sqrt(sum);
}

template <typename T, int D>
std::string ToString(Vec<T, D> value)
{
    std::ostringstream oss;
    for (int i = 0; i < D; i++)
    {
        oss << i << ' ';
    }
    return oss.str();
}

#endif /* defined(__Feynman_Simulator__vector__) */
