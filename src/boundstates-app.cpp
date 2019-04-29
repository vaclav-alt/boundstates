#include "boundstates-app.hpp"

namespace boundstates {

int BoundstatesApplication::Exec(int argc, const char **argv) {
    if (ParseCommandLineArguments(argc, argv)) {
		Initialize();
		InitLogging();
        Run();
	}
	return 0;
}

void BoundstatesApplication::Initialize() {
	CheckSettings();
}

int BoundstatesApplication::Run() {
    using namespace comptools;
    std::cout.precision(18);

    int basisSize = 100;
    double mu = 14682.6;
    // double mu = 1;

    auto V = input::ReadRRDFunctionFromFile(settings_.potentialFile, 1, 2);
    rmin_ = 0;
    rmax_ = 20;

    basis_ = basis::FourierBasis(rmax_, rmin_, basisSize);

    Eigen::MatrixXd mat(basisSize, basisSize);
    auto Vspline = interpol::NaturalSplineInterpol(V);

    for (int i = 0; i < basisSize; ++i) {
        for (int j = 0; j < basisSize; ++j) {
            mat(i,j) = MatrixElement(i,j, Vspline);
        }
        auto f = basis_.Freq(i);
        mat(i,i) += f * f / (2 * mu); 
    }

    // std::cout << mat << std::endl;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(mat);
    std::cout << std::fixed << eigensolver.eigenvalues() << std::endl;
}

double BoundstatesApplication::MatrixElement(int i, int j, std::function<double(double)> f) {
    using namespace comptools;

    auto romberg = integrator::Romberg<double>(rmin_, rmax_, [i, j,this, &f](double x) {
                return this->basis_(i,x) * f(x) * this->basis_(j,x);
            });
    
    return romberg.Integrate(10);

}

/* --- PROGRAM OPTIONS  --- */

boost::program_options::options_description BoundstatesApplication::DefineCommandLineArguments() {
    namespace po = boost::program_options;

	po::options_description desc("Allowed options");

	desc.add_options()
		("help", "produce help message")
        ("potential-file,f", po::value<PathType>(&settings_.potentialFile)->required(), "potential file")
		("output-path,o", po::value<PathType>(&settings_.outputPath)->default_value("./output/"), "output path")
		("wavefunctions,w", po::value<bool>(&settings_.calcwf)->default_value(0), "calculate wavefunctions")
		("logfile,l", po::value<PathType>(&settings_.logfile)->default_value(PathType("log")), "logfile")
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
			std::cout << "Usage: options_description [options]\n";
			std::cout << desc;
            return false;
        }

		po::notify(vm);

    }
    catch(std::exception& e) {
		std::cout << e.what() << std::endl;
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

	if ( settings_.logfile.parent_path().string() == "") {
		settings_.logfile = settings_.outputPath / settings_.logfile;
	}
}

/* --- LOGGING --- */

void BoundstatesApplication::InitLogging() {
	using namespace boost::log;
	register_simple_formatter_factory<trivial::severity_level, char>("Severity");
	add_file_log(
		keywords::file_name = settings_.logfile,
		keywords::format = "[%TimeStamp%] [%ThreadID%] [%Severity%] %Message%"
	);

	core::get()->set_filter (
		trivial::severity >= trivial::trace
	);

	add_common_attributes();
}

} // namespace boundstates
