#include <math.h>
#include "Vector3f.h"
#include "parser.h"

Vector3f Vector3f::normalize() {
   return (*this) /= this->length();
}

Vector3f Vector3f::cross(Vector3f const & v) const {
   return Vector3f(y*v.z - v.y*z, v.x*z - x*v.z, x*v.y - v.x*y);
}

double Vector3f::dot(Vector3f const & v) const {
   return x*v.x + y*v.y + z*v.z;
}

double Vector3f::length() const {
   return sqrtf(x*x + y*y + z*z);
}

Vector3f Vector3f::operator + (Vector3f const & v) const {
   return Vector3f(x+v.x, y+v.y, z+v.z);
}

Vector3f & Vector3f::operator += (Vector3f const & v) {
   x += v.x;
   y += v.y;
   z += v.z;

   return * this;
}

Vector3f Vector3f::operator - (Vector3f const & v) const {
   return Vector3f(x-v.x, y-v.y, z-v.z);
}

Vector3f & Vector3f::operator -= (Vector3f const & v) {
   x -= v.x;
   y -= v.y;
   z -= v.z;

   return * this;
}

Vector3f Vector3f::operator * (Vector3f const & v) const {
   return Vector3f(x*v.x, y*v.y, z*v.z);
}

Vector3f & Vector3f::operator *= (Vector3f const & v) {
   x *= v.x;
   y *= v.y;
   z *= v.z;

   return * this;
}

Vector3f Vector3f::operator / (Vector3f const & v) const {
   return Vector3f(x/v.x, y/v.y, z/v.z);
}

Vector3f & Vector3f::operator /= (Vector3f const & v) {
   x /= v.x;
   y /= v.y;
   z /= v.z;

   return * this;
}

Vector3f Vector3f::operator * (double const s) const {
   return Vector3f(x*s, y*s, z*s);
}

Vector3f & Vector3f::operator *= (double const s) {
   x *= s;
   y *= s;
   z *= s;

   return * this;
}

Vector3f Vector3f::operator / (double const s) const {
   return Vector3f(x/s, y/s, z/s);
}

Vector3f & Vector3f::operator /= (double const s) {
   x /= s;
   y /= s;
   z /= s;

   return * this;
}
