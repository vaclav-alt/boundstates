#include <cstdio>
#include <functional>

#include <comptools/function.hpp>
#include <comptools/interpolation.hpp>

#include "schrodinger-box.hpp"


class LcpApp {
public:
	struct Settings {
		double a = 1.7;
		double b = 10.0;
		double mu = 14582.6;
		size_t basisSize = 200;
		std::string V0file = "V0.dat";
	};

	void Exec(int argc, const char* argv[]) { Calculate(); };
	void Calculate() {
		auto VibrationalState = GetVibrationalState(0);
		std::FILE * fp = std::fopen("koule.dat", "w");
		for (auto p : VibrationalState) {
			fprintf(fp, "%.10f %.10f %.10f\n", p.x, p.y.real(), p.y.imag());
		}
		std::fclose(fp);

		SolveLcp();
		ProcessOutput();
	}
	comptools::Function<double, std::complex<double>> GetVibrationalState(size_t nu) {
		auto V = comptools::LoadFunctionFromFile(settings_.V0file);
		auto Vinterpolated = comptools::interpol::NaturalSplineInterpol(V);
		
		Parameters p;
		p.a = setting_.a;
		p.b = settings_.b;
		p.mu = settings_.mu;
		p.basisSize = settings_.basisSize;
		p.V = std::bind(Vinterpolated, std::placeholders::_1);
		
		SchrodingerBox box(p);
		printf("E[0] = %.10f\n", box.Energies()(0));
		
		comptools::Function<double, std::complex<double>> vibstate;
		for (size_t i = 0; i < 1001; ++i) {
			double x = 1.7 + i *(10.0 - 1.7) / 1000;
			comptools::Function<double, std::complex<double>>::Point p = { x, box.wf(0, x) };
			vibstate.push_back(p);
		}

		return vibstate;
	}

	void SolveLcp() {}
	void ProcessOutput() {}
private:
	Settings settings_;
};
