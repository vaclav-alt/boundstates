#pragma once

#include <stdexcept>
#include <iostream>
#include <cstdio>
#include <functional>

#include <comptools/interpolation.hpp>
#include <comptools/integration.hpp>
#include <comptools/function.hpp>
#include <comptools/grid.hpp>

#include <Eigen/Dense>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

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
	using MatrixType = Eigen::MatrixXd;
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
	void Diagonalize(MatrixType &);
	void EvaluatePotential();
	void TransformBack();

	double HamMatPontentialDiag(size_t i);
	double HamMatPontentialOffDiag(size_t i, size_t j);

    double MatrixElement(int i, int j, std::function<double(double)> f);

	void Initialize();
	void CheckSettings();

    bool ParseCommandLineArguments(int argc, const char **argv);
    int Run();
};


} // namespace fnscat
