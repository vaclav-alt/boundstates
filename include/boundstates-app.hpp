#pragma once

#include <cstdio>
#include "schrodinger-box.hpp"

// #include <boost/filesystem.hpp>

// #include <boost/program_options.hpp>
// #include <boost/property_tree/ptree.hpp>

namespace boundstates {

// typedef boost::filesystem::path PathType;

// namespace bpt = boost::property_tree;

// struct Settings {
// 	PathType outputPath;
// 	PathType potentialFile;
//     bool calcwf;
// };

class BoundstatesApplication {
public:
    BoundstatesApplication() {}
    
	int Exec(int argc, const char **argv);

private:
	// Settings settings_;
    
	// boost::program_options::options_description DefineCommandLineArguments();
	//void HamMatInitialize();

    // bool ParseCommandLineArguments(int argc, const char **argv);
};


} // namespace boudnstates
