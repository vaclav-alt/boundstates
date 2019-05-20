#include "boundstates-app.hpp"

int BoundstatesApplication::Exec(int argc, const char **argv) {
    if (ParseCommandLineArguments(argc, argv)) {
		using namespace std::placeholders;
		Parameters p;
		auto V = comptools::LoadFunctionFromFile(settings_.potentialFile.string());
		auto Vspline = comptools::interpol::NaturalSplineInterpol(V);

		p.basisSize = settings_.N;
		// p.mu = 14582.6;
		p.mu = settings_.mu;
		if (settings_.a == 0 && settings_.b == 0) {
			p.a = V.x(0);
			p.b = V.x(V.Size()-1);
			settings_.a = p.a;
			settings_.b = p.b;
		}
		else {
			p.a = settings_.a;
			p.b = settings_.b;
		}
		p.V = std::bind(Vspline, std::placeholders::_1);
		
		SchrodingerBox box(p);
		pbox = &box;
		// box.Run();

		SaveEigenenergiesForGnuplot(box.Energies());
		SaveAllWavefunctions();
		SaveInterpolatedPotential(Vspline);
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
        (",N",
		 	po::value<size_t>(&settings_.N)->default_value(100),
			"number of basis functions")
        ("mass,m",
		 	po::value<double>(&settings_.mu)->default_value(1.0),
			"mass")
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

void BoundstatesApplication::SaveEigenenergiesForGnuplot(SchrodingerBox::VectorType & energies) {
	std::FILE * fout = std::fopen("energies.dat", "w");
	size_t i = 1;
	fprintf(fout, "array E[%lu]\n", energies.size());
	for (auto & x : energies) 
		fprintf(fout, "E[%lu]=%12.8f\n", i++, x);
	fclose(fout);
}

void BoundstatesApplication::SaveWavefunction(size_t i) const {
	char filename[10];
	sprintf(filename, "wf%03lu.dat", i);
	std::FILE * fout = std::fopen(filename, "w");
	double h = (settings_.b - settings_.a) / 1000;
	for (int j = 0; j < 1000; ++j) {
		double x = settings_.a + j * h;
		double y = pbox->wf(i, x);
		fprintf(fout, "%12.8f %12.8f\n", x, y);
	}
	fclose(fout);
}

void BoundstatesApplication::SaveAllWavefunctions() const {
	for (size_t i = 0; i < settings_.N; ++i)
		SaveWavefunction(i);
}

void BoundstatesApplication::SaveInterpolatedPotential(const std::function<double(double)> & V) const {
	char const * filename = "V0_interpolated.dat";
	std::FILE * fout = std::fopen(filename, "w");
	double h = (settings_.b - settings_.a) / 10000;
	for (int j = 0; j < 10000; ++j) {
		double x = settings_.a + j * h;
		double y = V(x);
		fprintf(fout, "%12.8f %12.8f\n", x, y);
	}
	fclose(fout);
}
