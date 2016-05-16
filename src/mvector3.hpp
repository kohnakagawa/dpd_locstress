#pragma once

#include <cmath>
#include <iostream>

#ifdef __CUDACC__
#define INLINE __host__ __device__  __forceinline__ 
#else
#define INLINE __attribute__((always_inline))
#endif

//T should be double or float
template<typename T>
struct vector3{
  T x,y,z;
  
  //constructors
  INLINE vector3<T>() : x(0.0), y(0.0), z(0.0) {}
  INLINE vector3<T>(const vector3<T> &rhs) : x(rhs.x), y(rhs.y), z(rhs.z) {}
  INLINE explicit vector3<T>(const T &r) : x(r), y(r), z(r){}
  INLINE explicit vector3<T>(const T *p){
    x = p[0]; y = p[1]; z = p[2];
  }
  INLINE vector3<T>(const T &x_,const T &y_,const T &z_) : x(x_), y(y_), z(z_){}
 
  //operators
  INLINE const vector3<T> operator+(const vector3<T>& obj) const {
    const vector3<T> ret(x+obj.x,y+obj.y,z+obj.z);
    return ret;
  }
  INLINE const vector3<T> operator-(const vector3<T>& obj) const {
    const vector3<T> ret(x-obj.x,y-obj.y,z-obj.z);
    return ret;
  }
  INLINE const vector3<T>& operator+=(const vector3<T>& obj){
    x += obj.x; y += obj.y; z += obj.z;
    return *this;
  }
  INLINE const vector3<T>& operator-=(const vector3<T>& obj){
    x -= obj.x; y -= obj.y; z -= obj.z;
    return *this;
  }
  INLINE const vector3<T> operator*(const T& cf) const {
    const vector3<T> ret(x*cf,y*cf,z*cf);
    return ret;
  }
  INLINE const T operator*(const vector3<T>& obj) const {
    const T ret = x*obj.x + y*obj.y + z*obj.z;
    return ret;
  }
  INLINE friend const vector3<T> operator*(const T& cf,const vector3<T>& obj) {
    const vector3<T> ret(obj.x*cf,obj.y*cf,obj.z*cf);
    return ret;
  }
  INLINE const vector3<T>& operator*=(const T& cf){
    x *= cf; y *= cf; z *= cf;
    return *this;
  }
  INLINE const vector3<T>& operator/=(const T& cf){
    x /= cf; y /= cf; z /= cf;
    return *this;
  }
  INLINE friend const vector3<T> operator/(const T& cf,const vector3<T>& obj) {
    const vector3<T> ret(cf/obj.x,cf/obj.y,cf/obj.z);
    return ret;
  }
  INLINE T &operator [](int i){
    return (&x)[i];
  }
  INLINE bool operator<(const vector3<T>& rhs)const{
    const bool ret = (x < rhs.x) && (y < rhs.y) && (z < rhs.z);
    return ret;
  }
  INLINE bool operator>(const vector3<T>& rhs)const{
    const bool ret = (x > rhs.x) && (y > rhs.y) && (z > rhs.z);
    return ret;
  }
  INLINE bool operator<=(const vector3<T>& rhs)const{
    const bool ret = (x <= rhs.x) && (y <= rhs.y) && (z <= rhs.z);
    return ret;
  }
  INLINE bool operator>=(const vector3<T>& rhs)const{
    const bool ret = (x >= rhs.x) && (y >= rhs.y) && (z >= rhs.z);
    return ret;
  }
  
  //inverse square root
  //for IEEE754
  INLINE float m_rsqrt(const float& r) const {
    const float rHalf = 0.5 * r;
    const int   tmp   = 0x5F3759DF - ( *(int*)&r >> 1 );
    float rRes        = *(float*)&tmp;
    rRes *= ( 1.5 - ( rHalf * rRes * rRes ) );
    rRes *= ( 1.5 - ( rHalf * rRes * rRes ) );
    return rRes;
  }
  //for IEEE754
  INLINE double m_rsqrt(const double& r) const {
    const double         rHalf = 0.5 * r;
    const long long int  tmp   = 0x5FE6EB50C7B537AAl - ( *(long long int*)&r >> 1);//initial guess
    double         rRes        = *(double*)&tmp;
    rRes *= ( 1.5 - ( rHalf * rRes * rRes ) );
    rRes *= ( 1.5 - ( rHalf * rRes * rRes ) );
    return rRes;
  }
  
  //ret norm
  INLINE const T invnorm2() const {
    const T dr2 = (*this)*(*this);
    return m_rsqrt(dr2);
  }
  INLINE const T norm2() const {
    const T dr2 = (*this)*(*this);
    return sqrt(dr2);
  }
  INLINE const T fastnorm2() const {
    const T dr2 = (*this)*(*this);
    return m_rsqrt(dr2) * dr2;
  }
  INLINE const T dist2() const {
    const T dr2 = (*this)*(*this);
    return dr2;
  }
  INLINE void clear() {
    x = y = z = 0.0;
  }
  
#ifndef __CUDACC__
  INLINE bool isfinite3() const {
    const bool ret = (std::isfinite(x)) && (std::isfinite(y)) && (std::isfinite(z));
    if(!ret) std::cout << x << " " << y << " " << z << std::endl;
    return ret;
  }
#endif

  friend std::ostream& operator<<(std::ostream& os, const vector3& dvec3){
    os << dvec3.x << " " << dvec3.y << " " << dvec3.z;
    return os;
  }
    
};

#ifndef __CUDACC__
typedef vector3<double> double3;
typedef vector3<float > float3;
#else
typedef vector3<double> dbl3;
typedef vector3<float > flt3;
#endif

#undef INLINE
