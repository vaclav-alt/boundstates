#pragma once

#include <cstdio>
#include <functional>

#include <comptools/array.hpp>
#include <comptools/function.hpp>
#include <comptools/grid.hpp>
#include <comptools/interpolation.hpp>

#include <lapacke.h>
#include <cblas.h>

struct Parameters {
	size_t basisSize = 50;
	double a = -10;
	double b = 10;
	double mu = 1.0;
	std::function<double(double)> V = [](double x) { return x*x/2; };
};

class SchrodingerBox {
public:
	using MatrixType = comptools::Array<double,2>;
	using VectorType = comptools::Array<double,1>;
    SchrodingerBox() {}
    
    void Run();
private:
	Parameters p;
    
	void HamMatPotential(MatrixType &);
	void HamMatKinetic(MatrixType &);

	void FillWithXMean(MatrixType &);
	VectorType Diagonalize(MatrixType &);
	void TransformBack(MatrixType &, VectorType &);
	void EvaluatePotential(VectorType &);

	double PontentialMatrixElement(size_t i, size_t j);

	void printMatrix(MatrixType & m);
};
