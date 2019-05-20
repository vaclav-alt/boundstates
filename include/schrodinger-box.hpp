#pragma once

#include <cstdio>
#include <functional>

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
	using MatrixType = comptools::Array<double,2>;
	using VectorType = comptools::Array<double,1>;
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

   	double wf(size_t i, double x); 
	VectorType & Energies() { return energies; }
private:
	Parameters p;
	MatrixType hmat;
	VectorType energies;
    
    void Run();
	void HamMatPotential(MatrixType &);
	void HamMatKinetic(MatrixType &);

	void FillWithXMean(MatrixType &);
	VectorType Diagonalize(MatrixType &);
	void TransformBack(MatrixType &, VectorType &);
	void EvaluatePotential(VectorType &);

	double PontentialMatrixElement(size_t i, size_t j);
	double basis(size_t i, double x) {
		double L = p.b - p.a;
		return sqrt(2 / L) * sin(M_PI * (i+1) * (x - p.a) / L);
	}

	void printMatrix(MatrixType & m);
};
