#include "boundstates-app.hpp"

namespace boundstates {

int BoundstatesApplication::Exec(int argc, const char **argv) {
    if (ParseCommandLineArguments(argc, argv)) {
		Initialize();
        Run();
	}
	return 0;
}

void BoundstatesApplication::Initialize() {
	CheckSettings();
}

int BoundstatesApplication::Run() {

	MatrixType hmat(p.basisSize, p.basisSize);
	HamMatPotential(hmat);
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(mat);
    // std::cout << std::fixed << eigensolver.eigenvalues() << std::endl;
	HamMatKinetic(hmat);
	Eigen::SelfAdjointEigenSolver<MatrixType> es;
	es.compute(hmat);
	std::cout << es.eigenvalues() << std::endl;
	return 0;
}

void BoundstatesApplication::HamMatPotential(MatrixType & m) {
	FillWithXMean(m);
	Diagonalize(m);
	// EvaluatePotential();
	// TransformBack();
}

void BoundstatesApplication::FillWithXMean(MatrixType & m) {
	for (size_t i = 0; i < m.rows(); ++i) {
		for (size_t j = 0; j < i; ++j) {
			if ((i+j+1)%2 == 0) {
				m(i, j) = HamMatPontentialOffDiag(i+1 , j+1);
			} else {
				m(i,j) = 0;
			}
		}
		m(i,i) = 0;
	}
}

double BoundstatesApplication::HamMatPontentialDiag(size_t i) {
	double pisq = M_PI*M_PI;
	double Vdiag = -p.a + p.b + 8 * (p.a + p.b) * i * i * pisq
					+ (p.a - p.b) * cos(4 * i * M_PI)
					- 4 * p.b * i * M_PI * sin(4 * i * M_PI);
	// Vdiag *= (p.b - p.a) / (32 * i * i * pisq);
	Vdiag *= 2.0 / (32 * i * i * pisq);
	return Vdiag;
}

double BoundstatesApplication::HamMatPontentialOffDiag(size_t i, size_t j) {
	size_t dif = i - j;
	size_t sum = i + j;
	double Vdiag = -8.0 * i * j * (p.b - p.a) / (M_PI*M_PI*sum*sum*dif*dif);
	return Vdiag;
}

void BoundstatesApplication::Diagonalize(MatrixType & m) {
	Eigen::SelfAdjointEigenSolver<MatrixType> es;
	es.compute(m);
	auto x = es.eigenvalues();
	MatrixType V = es.eigenvectors();
	for (size_t i = 0; i < x.size(); ++i) {
		x[i] = p.V(x[i]);
	}
	MatrixType D = x.asDiagonal();
	m = V * D * V.transpose();
}

void BoundstatesApplication::HamMatKinetic(MatrixType & m) {
	double f1 = M_PI / (p.b - p.a);
	for (size_t i = 0; i < m.rows(); ++i) {
		double f = (i+1) * f1;
		m(i, i) += f * f / 2 / p.mu;
	}
}

double BoundstatesApplication::MatrixElement(int i, int j, std::function<double(double)> f) {
    using namespace comptools;

    return Integrate<double>::Romberg(rmin_, rmax_, f, 10);

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
}

} // namespace boundstates
