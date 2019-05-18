#pragma once

#include <boost/filesystem.hpp>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

#include <cstdio>
#include <iostream>
#include <functional>
#include <stdexcept>

#include "schrodinger-box.hpp"
#include "comptools/function.hpp"
#include "comptools/interpolation.hpp"

typedef boost::filesystem::path PathType;

namespace bpt = boost::property_tree;

struct Settings {
	double a;
	double b;
	PathType outputPath;
	PathType potentialFile;
    bool calcwf;
};

class BoundstatesApplication {
public:
    BoundstatesApplication() {}
    
	int Exec(int argc, const char **argv);

private:
	Settings settings_;
    
	void CheckSettings();
	boost::program_options::options_description DefineCommandLineArguments();
    bool ParseCommandLineArguments(int argc, const char **argv);
};
