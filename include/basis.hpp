#pragma once

#include <cmath>

namespace basis {

class BasisBase {
public:
    BasisBase(int size) :
        size_(size)
    {}

    virtual double operator()(int i, double x) {}
    virtual double d(int i, double x) {}
    virtual double d2(int i, double x) {}

    int Size() const {
        return size_;
    }
private:
    int size_;
};


class SineBasis : public BasisBase {
public:

    SineBasis() :
        SineBasis(1,0,0)
    {}

    SineBasis(double Rmax, double Rmin, int size) :
        L_(Rmax - Rmin),
        Rmin_(Rmin),
        BasisBase(size)
    {}

    double operator()(int i, double R) {
        double x = R - Rmin_;
        auto f = Freq(i);
        return std::sin(f * x) * Norm();
    }

    double Freq(int i) {
        return (i+1) * M_PI / L_;
    }


private:
    double Norm() {
        return std::sqrt(2 / L_);
    }


    double L_;
    double Rmin_;
};

class FourierBasis : public BasisBase {
public:

    FourierBasis() :
        FourierBasis(1,0,0)
    {}

    FourierBasis(double Rmax, double Rmin, int size) :
        L_(Rmax - Rmin),
        Rmin_(Rmin),
        BasisBase(size)
    {}

    double operator()(int i, double R) {
        double x = R - Rmin_;
        auto f = Freq(i);

        if ( (i % 2) == 0){
            return std::cos(f * x) * Norm();
        }
        else {
            return std::sin(f * x) * Norm();
        }
    }

    double d(int i, double x) {
        auto f = Freq(i);
        if ( (i % 2) == 0){
            return f * std::sin(f * x) * Norm();
        }
        else {
            return f * std::cos(f * x) * Norm();
        }
    }

    double d2(int i, double x) {
        auto f = Freq(i);
        return - f * f * (*this)(i,x);
    }

    double Freq(int i) {
        int j = (i+1)/2;
        return j * M_PI / L_;
    }


private:
    double Norm() {
        return std::sqrt(2 / L_);
    }


    double L_;
    double Rmin_;
};

}  // namespace basis
