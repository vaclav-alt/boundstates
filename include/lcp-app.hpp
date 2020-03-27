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
		double c = 8.5;
		double eta = M_PI/8;
		double ea = 0.053695788;
		double mu = 14582.6;
		size_t NGridPoints = 8300;
		size_t basisSize = 200;
		size_t nu = 0;
		std::string V0file = "V0.dat";
		std::string Vlocfile = "Vloc.dat";
	};
	struct RuntimeResults {
		double Enu = 0.0;
		double VresAs = 0.0;

	};

	void Exec(int argc, const char* argv[]) { Calculate(); };
	TridiagonalSystem GenerateTridiagonalSystem();
	void UpdateEnergyTerm(double E, TridiagonalSystem &);
	void Calculate() {

		xgrid_ = comptools::Grid::FromMinMaxN(settings_.a,
														settings_.b,
														settings_.NGridPoints);
		auto Egrid = comptools::Grid::FromMinMaxN(0.0, 0.50, 5000);
		ecsgrid_ = comptools::EcsGrid::FromMinMaxN(settings_.a,
														settings_.b,
														settings_.NGridPoints,
														settings_.c,
														settings_.eta);

		PrecalculateFunctions();
		auto system = GenerateTridiagonalSystem();

		FILE * fp = std::fopen("csda.dat", "w");

		size_t i = 0;
		for (auto E : Egrid) {
			auto eqsystem = system;
			UpdateEnergyTerm(E, eqsystem);
			Thomas(eqsystem);

			fprintf(fp, "%.10f %.10f\n", E, CrossSection(eqsystem, E, ecsgrid_.ic()-1));
			if ((i%100)==0) {
				char filename[15];
				sprintf(filename, "LCPwf%04lu.dat", i);
				FILE * fp2 = std::fopen(filename, "w");
				for (size_t i = 0; i < ecsgrid_.size(); ++i) {
					fprintf(fp2, "%.10f %.10f %.10f\n", xgrid_[i], eqsystem.b[i].real(), eqsystem.b[i].imag());
				}
				std::fclose(fp2);
			}
			++i;
		}
		std::fclose(fp);

		SaveFunctions();
	}

	void GetVibrationalState(size_t);
	double CrossSection(TridiagonalSystem & s, double E, size_t i) {
		double phisq= std::abs(s.b[i]);
		double K = sqrt(2 * settings_.mu * (E - rr_.VresAs));
		double kvsq = 2 * (E - rr_.Enu);
		return 2 * M_PI * M_PI / kvsq * K / settings_.mu * phisq;
	}

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
	RuntimeResults rr_;
	Functions functions_;
	comptools::Grid xgrid_;
	comptools::EcsGrid ecsgrid_;
};
