#include "boundstates-app.hpp"

namespace boundstates {

int BoundstatesApplication::Exec(int argc, const char **argv) {
//    if (ParseCommandLineArguments(argc, argv)) {
	SchrodingerBox box;
	box.Run();
//	}
	return 0;
}


/* --- PROGRAM OPTIONS  --- */
// 
// boost::program_options::options_description BoundstatesApplication::DefineCommandLineArguments() {
//     namespace po = boost::program_options;
// 
// 	po::options_description desc("Allowed options");
// 
// 	desc.add_options()
// 		("help", "produce help message")
//         ("potential-file,f", po::value<PathType>(&settings_.potentialFile)->required(), "potential file")
// 		("output-path,o", po::value<PathType>(&settings_.outputPath)->default_value("./output/"), "output path")
// 		("wavefunctions,w", po::value<bool>(&settings_.calcwf)->default_value(0), "calculate wavefunctions")
// 	;
// 
// 	return desc;
// }
// 
// bool BoundstatesApplication::ParseCommandLineArguments(int argc, const char **argv) {
//     namespace po = boost::program_options;
// 
//     try {
// 		auto desc = DefineCommandLineArguments();
// 
//         po::variables_map vm;
//         po::store(po::parse_command_line(argc, argv, desc), vm);
// 
//         if (vm.count("help")) {
// 			std::cout << "Usage: options_description [options]\n";
// 			std::cout << desc;
//             return false;
//         }
// 
// 		po::notify(vm);
// 
//     }
//     catch(std::exception& e) {
// 		std::cout << e.what() << std::endl;
//         return false;
//     }
//     return true;
// }
// 
// void BoundstatesApplication::CheckSettings() {
// 	namespace bf = boost::filesystem;
// 
// 	if (!bf::exists(settings_.outputPath)) {
// 		bf::create_directory(settings_.outputPath);
// 	}
// 	else if (!bf::is_directory(settings_.outputPath)) {
// 		throw std::runtime_error("output-path is not a directory: " + settings_.outputPath.string());
// 	}
// }
// 
} // namespace boundstates
