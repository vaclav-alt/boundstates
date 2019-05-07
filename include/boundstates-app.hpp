#pragma once

#include <stdexcept>
#include <iostream>
#include <functional>

#include <comptools/interpolation.hpp>
#include <comptools/integration.hpp>
#include <comptools/function.hpp>
#include "basis.hpp"
#include <comptools/grid.hpp>

#include <Eigen/Dense>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <boost/log/core.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>


namespace boundstates {

typedef boost::filesystem::path PathType;

namespace bpt = boost::property_tree;

struct Settings {
	PathType outputPath;
	PathType logfile;
	PathType potentialFile;
    bool calcwf;
};

class BoundstatesApplication {
public:
    BoundstatesApplication() {}
    
	int Exec(int argc, const char **argv);

private:
	Settings settings_;
    basis::SineBasis basis_;
    double rmin_ = 0;
    double rmax_ = 1;

    
	boost::program_options::options_description DefineCommandLineArguments();

    double MatrixElement(int i, int j, std::function<double(double)> f);
    void SetUpBasis();

    void Eigenvalues();
    void Eigenproblem();

	void Initialize();
	void CheckSettings();
	void InitLogging();
    bool ParseCommandLineArguments(int argc, const char **argv);
    int Run();
};


} // namespace fnscat
