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
	HamMatKinetic(hmat);
	auto D = Diagonalize(hmat);
	for (auto & x : D) 
		printf("%12.8f\n", x);
	return 0;
}

void BoundstatesApplication::HamMatPotential(MatrixType & m) {
	FillWithXMean(m);
	auto D = Diagonalize(m);
	EvaluatePotential(D);

	TransformBack(m, D);
}

void BoundstatesApplication::printMatrix(MatrixType & m) {
	for (size_t i = 0; i < m.rows(); ++i) {
		for (size_t j = 0; j < m.cols(); ++j) {
			printf("%12.8f ", m(i,j));
		}
		printf("\n");
	}
}

void BoundstatesApplication::FillWithXMean(MatrixType & m) {
	for (size_t i = 0; i < m.rows(); ++i) {
		for (size_t j = 0; j < i; ++j) {
			if ((i+j+1)%2 == 0) {
				m(i, j) = PontentialMatrixElement(i+1 , j+1);
			} else {
				m(i,j) = 0;
			}
		}
		m(i,i) = 0;
	}
}

double BoundstatesApplication::PontentialMatrixElement(size_t i, size_t j) {
	size_t dif = i - j;
	size_t sum = i + j;
	double Vdiag = -8.0 * i * j * (p.b - p.a) / (M_PI*M_PI*sum*sum*dif*dif);
	return Vdiag;
}

BoundstatesApplication::VectorType BoundstatesApplication::Diagonalize(MatrixType & A) {
	char JOBZ = 'V';
	char UPLO = 'U';
	int N = A.rows();
	VectorType W(N);
	int LWORK=3*N-1;
	std::vector<double> WORK(LWORK);
	int INFO;

	dsyev_(&JOBZ,
			&UPLO,
			&N,
			A.data(),
			&N,
			W.data(),
			WORK.data(),
			&LWORK,
			&INFO);

	return W;
}

void BoundstatesApplication::EvaluatePotential(VectorType & v) {
	for (auto & x : v)
		x = p.V(x);
}

void BoundstatesApplication::TransformBack(MatrixType & m, VectorType & D) {
	auto A = m;
	for (size_t i = 0; i < m.rows(); ++i) {
		for (size_t j = 0; j < m.cols(); ++j) {
			double temp = 0;
			for (size_t k = 0; k < m.cols(); ++k) {
				temp += A(k,i) * A(k,j) * D(k);
			}
			m(i,j) = temp;
		}
	}
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
