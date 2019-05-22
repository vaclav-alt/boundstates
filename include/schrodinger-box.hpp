#pragma once

#include <cstdio>
#include <functional>
#include <complex.h>

#define lapack_complex_double std::complex<double>

#include <comptools/array.hpp>
#include <comptools/function.hpp>
#include <comptools/grid.hpp>

#include <lapacke.h>
#include <cblas.h>


struct Parameters {
	size_t basisSize = 100;
	double a = 1.7;
	double b = 10;
	double mu = 1.0;
	std::function<double(double)> V = [](double x) { return x*x/2; };
};

class SchrodingerBox {
public:
	using NumberType = lapack_complex_double;
	using MatrixType = comptools::Array<NumberType,2>;
	using VectorType = comptools::Array<NumberType,1>;
	using RealVectorType = comptools::Array<double,1>;
    SchrodingerBox() :
		SchrodingerBox(Parameters())
	{}
    SchrodingerBox(Parameters p) :
		p(p),
		hmat(p.basisSize, p.basisSize),
		energies(p.basisSize)
	{
		Run();
	}

   	NumberType wf(size_t i, double x); 
	RealVectorType & Energies() { return energies; }
private:
	Parameters p;
	MatrixType hmat;
	RealVectorType energies;
    
    void Run();
	void HamMatPotential(MatrixType &);
	void HamMatKinetic(MatrixType &);

	void FillWithXMean(MatrixType &);
	RealVectorType Diagonalize(MatrixType &);
	void TransformBack(MatrixType &, RealVectorType &);
	void EvaluatePotential(RealVectorType &);

	NumberType PontentialMatrixElement(size_t i, size_t j);
	NumberType basis(size_t i, double x) {
		double L = p.b - p.a;
		double f = M_PI * (i+1) / L;
		return sqrt(2 / L) * sin(f * (x - p.a));
	}

	void printMatrix(MatrixType & m);
};
