#pragma once

#include <cmath>
#include <complex>

#include <comptools/math.hpp>

class DomckeModel {
public:
	struct Parameters {
		double al = 1.25;
		double A = 1.757;
		double B = 1.667;
		double C = 0.98;
		double mu = 1.0;
		double a = 1.96;
		double D0 = 5.0;
		double D1 = 3.2;
		double R0 = 0.0;
		double t = 0.2;
		double Qd = 1.2;
		double bisa = -5;
		double bisb = 10;
	};
	
	DomckeModel(Parameters p) :
		p(p)
	{}

	double f(double E) {
		return sqrt(0.5 * p.A / M_PI) * std::pow(E / p.B, p.al) * exp(-0.5 * E / p.B);
	}

	double g(double R) {
		auto x = R - p.R0 + 0.5 / sqrt(p.C);
		return exp(-p.C * x*x);
	}

	double Gamma(double E) {
		auto fE = f(E);
		return 2 * M_PI * fE * fE;
	}

	double Vde(double E, double R) {
		return f(E) * g(R);
	}

	double Gamma(double E, double R) {
		return 2 * M_PI * Vde(E, R);
	}

	double Delta(double E) {
		using comptools::math::DeltaIntegral;
		return 0.5 * p.A * std::pow(1.0/p.B, 2.0*p.al) / M_PI * DeltaIntegral(E, 2.0*p.al, 1.0/p.B);
	}

	double Delta(double E, double R) {
		auto gR = g(R);
		return gR * gR * Delta(E);
	}
	
	double Vres(double R) {
		return V0(R) + Eres(R);
	}

	double V0(double R) {
		auto tmp = exp(-p.a * (R - p.R0)) - 1;
		return p.D0 * tmp * tmp;
	}

	double Vd(double R) {
		auto tmp = exp(-p.a * (R - p.R0));
		return p.D1 * tmp * (tmp - 2 * p.t) + p.Qd;
	}
	
	std::complex<double> Vloc(double R) {
		using namespace std::complex_literals;
		auto eres = Eres(R);
		double gamma = (eres >= 0) ? Gamma(eres, R) : 0.0;
		return V0(R) + eres + gamma * 1i;
	}

	double EresEquation(double E, double R) {
		return Delta(E, R) + Vd(R) - V0(R) - E;
	}

private:
	double Eres(double R) {
		auto eres = comptools::math::Bisection([&R,this](double E) {
					return EresEquation(E,R);
				}, p.bisa, p.bisb);
		return eres;
	}

	Parameters p;
};
