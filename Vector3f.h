#ifndef __VECTOR3f_H__
#define __VECTOR3f_H__

class Vector3f {
public:

   double x, y, z;

   Vector3f() : x(0), y(0), z(0) {}

   Vector3f(double in) : x(in), y(in), z(in) {}

   Vector3f(double in_x, double in_y, double in_z) : x(in_x), y(in_y), z(in_z) {}

   Vector3f normalize();

   Vector3f cross(Vector3f const & v) const;

   double dot(Vector3f const & v) const;

   double length() const;

   Vector3f operator + (Vector3f const & v) const;

   Vector3f & operator += (Vector3f const & v);

   Vector3f operator - (Vector3f const & v) const;

   Vector3f & operator -= (Vector3f const & v);

   Vector3f operator * (Vector3f const & v) const;

   Vector3f & operator *= (Vector3f const & v);

   Vector3f operator / (Vector3f const & v) const;

   Vector3f & operator /= (Vector3f const & v);

   Vector3f operator * (double const s) const;

   Vector3f & operator *= (double const s);

   Vector3f operator / (double const s) const;

   Vector3f & operator /= (double const s);
};

#endif
