#include <cstdio>
#include <complex>
#include <functional>
#include <vector>

#include <comptools/function.hpp>
#include <comptools/ecs_grid.hpp>
#include <comptools/interpolation.hpp>

#include "schrodinger-box.hpp"
#include "csv.h"
#include "complex_thomas.hpp"

using RealFunctionType = comptools::Function<double, double>;
using ComplexFunctionType = comptools::Function<double, std::complex<double>>;

class LcpApp {
public:
	struct Functions {
		RealFunctionType vres;
		RealFunctionType gamma;
		ComplexFunctionType chi;
	};

	struct Settings {
		double a = 1.7;
		double b = 10.0;
		double c = 8;
		double eta = M_PI/8;
		double ea = -0.05379;
		double mu = 14582.6;
		double E = -0.05;
		size_t NGridPoints = 83000;
		size_t basisSize = 200;
		std::string V0file = "V0.dat";
		std::string Vlocfile = "Vloc.dat";
	};

	LcpApp(Settings s) :
		settings_(s)
	{}

	LcpApp()
	{}

	void Exec(int argc, const char* argv[]) { Calculate(); };
	TridiagonalSystem GenerateTridiagonalSystem();
	void UpdateEnergyTerm(double E, TridiagonalSystem &);
	void Calculate() {

		xgrid_ = comptools::Grid::FromMinMaxN(settings_.a,
														settings_.b,
														settings_.NGridPoints);
		ecsgrid_ = comptools::EcsGrid::FromMinMaxN(settings_.a,
														settings_.b,
														settings_.NGridPoints,
														settings_.c,
														settings_.eta);

		PrecalculateFunctions();
		auto system = GenerateTridiagonalSystem();
		UpdateEnergyTerm(settings_.E, system);
		Thomas(system);

		FILE * fp = std::fopen("wavefunction.dat", "w");
		for (size_t i = 0; i < ecsgrid_.size(); ++i) {
			fprintf(fp, "%.10f %.10f %.10f\n", xgrid_[i], system.b[i].real(), system.b[i].imag());
		}
		std::fclose(fp);

		SaveFunctions();
	}
	void GetVibrationalState(size_t);

	void PrecalculateFunctions();
	void SaveGridFtorToFile(comptools::Grid, const std::function<double(double)> &, std::string) const;
	void SaveFunctionToFile(ComplexFunctionType &, std::string) const;
	void SaveFunctionToFile(RealFunctionType &, std::string) const;
	void SaveFunctions() {
		SaveFunctionToFile(functions_.vres, "vres.dat");
		SaveFunctionToFile(functions_.gamma, "gamma.dat");
		SaveFunctionToFile(functions_.chi, "chi.dat");
	}
private:

	Settings settings_;
	Functions functions_;
	comptools::Grid xgrid_;
	comptools::EcsGrid ecsgrid_;
};
