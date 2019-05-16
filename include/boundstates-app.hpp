#pragma once

#include <stdexcept>
#include <iostream>
#include <cstdio>
#include <functional>

#include <comptools/array.hpp>
#include <comptools/interpolation.hpp>
#include <comptools/integration.hpp>
#include <comptools/function.hpp>
#include <comptools/grid.hpp>

#include <lapacke.h>
#include <cblas.h>

#include <boost/filesystem.hpp>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace boundstates {

typedef boost::filesystem::path PathType;

namespace bpt = boost::property_tree;

struct Settings {
	PathType outputPath;
	PathType potentialFile;
    bool calcwf;
};

struct Parameters {
	size_t basisSize = 50;
	double a = -10;
	double b = 10;
	double mu = 1.0;
	std::function<double(double)> V = [](double x) { return x*x/2; };
};

class BoundstatesApplication {
public:
	using MatrixType = comptools::Array<double,2>;
	using VectorType = comptools::Array<double,1>;
    BoundstatesApplication() {}
    
	int Exec(int argc, const char **argv);

private:
	Parameters p;
	Settings settings_;
    double rmin_ = 0;
    double rmax_ = 1;
    
	boost::program_options::options_description DefineCommandLineArguments();
	//void HamMatInitialize();
	void HamMatPotential(MatrixType &);
	void HamMatKinetic(MatrixType &);

	void FillWithXMean(MatrixType &);
	VectorType Diagonalize(MatrixType &);
	void TransformBack(MatrixType &, VectorType &);
	void EvaluatePotential(VectorType &);
	void TransformBack();

	double PontentialMatrixElement(size_t i, size_t j);

    double MatrixElement(int i, int j, std::function<double(double)> f);

	void Initialize();
	void CheckSettings();
	void printMatrix(MatrixType & m);

    bool ParseCommandLineArguments(int argc, const char **argv);
    int Run();
};


} // namespace fnscat
