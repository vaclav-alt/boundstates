#include "boundstates-app.hpp"

int BoundstatesApplication::Exec(int argc, const char **argv) {
    if (ParseCommandLineArguments(argc, argv)) {
		using namespace std::placeholders;
		Parameters p;
		auto V = comptools::LoadFunctionFromFile(settings_.potentialFile.string());
		auto Vspline = comptools::interpol::NaturalSplineInterpol(V);

		p.basisSize = 200;
		p.mu = 14582.6;
		if (settings_.a == 0 && settings_.b == 0) {
			p.a = V.x(0);
			p.b = V.x(V.Size()-1);
		}
		else {
			p.a = settings_.a;
			p.b = settings_.b;
		}
		p.V = std::bind(Vspline, std::placeholders::_1);
		
		SchrodingerBox box(p);
		box.Run();

		std::FILE * fout = std::fopen("Vspline.dat", "w");
		for (int i = 1; i < 1000; ++i) {
			double x = -5 + i * 0.01;
			double y = Vspline(x);
			fprintf(fout, "%12.8f %12.8f\n", x, y);
		}
		fclose(fout);
	}
	return 0;
}

/* --- PROGRAM OPTIONS  --- */

boost::program_options::options_description BoundstatesApplication::DefineCommandLineArguments() {
    namespace po = boost::program_options;

	po::options_description desc("Allowed options");

	desc.add_options()
		("help,h",
		 	"produce help message"
			)
        ("potential-file,f",
		 	po::value<PathType>(&settings_.potentialFile)->required(),
			"potential file")
        (",a",
		 	po::value<double>(&settings_.a)->default_value(0),
			"left boundary")
        (",b",
		 	po::value<double>(&settings_.b)->default_value(0),
			"right boundary")
		("output-path,o",
		 	po::value<PathType>(&settings_.outputPath)->default_value("./output/"),
			"output path")
		("wavefunctions,w",
		 	po::value<bool>(&settings_.calcwf)->default_value(0),
			"calculate wavefunctions")
	;

	return desc;
}

bool BoundstatesApplication::ParseCommandLineArguments(int argc, const char **argv) {
    namespace po = boost::program_options;

    try {
		auto desc = DefineCommandLineArguments();

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
			printf("Usage: options_description [options]\n");
			std::cout << desc << std::endl;
            return false;
        }

		po::notify(vm);

    }
    catch(std::exception& e) {
		printf("%s\n", e.what());
        return false;
    }
    return true;
}

void BoundstatesApplication::CheckSettings() {
	namespace bf = boost::filesystem;

	if (!bf::exists(settings_.outputPath)) {
		bf::create_directory(settings_.outputPath);
	}
	else if (!bf::is_directory(settings_.outputPath)) {
		throw std::runtime_error("output-path is not a directory: " + settings_.outputPath.string());
	}
}
