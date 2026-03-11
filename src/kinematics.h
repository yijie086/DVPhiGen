#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "event.h"
#include <string>

struct Vec3 {
    double x, y, z;

    Vec3();
    Vec3(double x_, double y_, double z_);

    double Mag2() const;
    double Mag() const;
    Vec3 Unit() const;

    Vec3 operator+(const Vec3& other) const;
    Vec3 operator-(const Vec3& other) const;
    Vec3 operator*(double a) const;
    Vec3 operator/(double a) const;
    double Dot(const Vec3& other) const;
    Vec3 Cross(const Vec3& other) const;
};

struct FourVector {
    double E, px, py, pz;

    FourVector();
    FourVector(double E_, double px_, double py_, double pz_);

    double M2() const;
    double M() const;
    Vec3 Vect() const;

    FourVector operator+(const FourVector& other) const;
    FourVector operator-(const FourVector& other) const;

    void Boost(const Vec3& beta);
};

struct FinalState {
    FourVector electron;
    FourVector proton;
    FourVector kp;
    FourVector km;

    void Print() const;
};

class KinematicsCalculator {
public:
    static FinalState BuildEvent(const DVphiEvent& evt);

private:
    static double Lambda(double x, double y, double z);
};

#endif